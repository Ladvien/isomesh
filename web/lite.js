// The front page's interactive demo: `isomesh_web.wasm` on one side, WebGL2 on
// the other, and nothing in between.
//
// # Why this is hand-written and has no dependencies
//
// The nine playable demos are Bevy builds and each is ~36 MB because it carries a
// renderer, an asset system and a window backend. This page is the other half of
// that argument: the module it drives is ~130 KB because `isomesh` is `no_std`
// with one dependency. Reaching for `wasm-bindgen` here would add 40-80 KB of
// glue to a 130 KB module and move the GL state machine into Rust, and reaching
// for a matrix library would add a package to a page whose whole claim is that it
// does not need one. So the ABI is `extern "C"` over linear memory and the six
// lines of matrix maths are below.
//
// # The one thing that will bite an editor of this file
//
// **Every typed-array view is re-created after every `iso_mesh` call.** Wasm
// memory growth detaches `memory.buffer`, so a `Float32Array` captured before a
// re-mesh reads a buffer that no longer exists -- silently, as zeros. The
// module's own docs say its pointers are valid only until the next `iso_mesh`,
// and `readMesh` below is the only place that reads them.
//
// # WebGL2, and there is no WebGL1 path
//
// `gl.drawElements(..., gl.UNSIGNED_INT, ...)` is core in WebGL2 and an extension
// in WebGL1, and a 49^3 Marching Tetrahedra mesh is past 65,535 vertices. A
// second renderer for a browser that cannot do this would be a second path to
// maintain and to get subtly wrong; if WebGL2 is missing this says so and stops.

const SAMPLES_DEFAULT = 33;

// The site's own `--panel`, so the canvas has no seam against the page.
const CLEAR = [0.063, 0.078, 0.11, 1.0];

// A single directional light, in view space, and a flat base colour. This is a
// demo of a surface, not of a shading model.
const LIGHT = [0.42, 0.78, 0.47];
const SURFACE = [0.62, 0.72, 0.85];
const WIRE = [0.98, 0.72, 0.31];

const VERTEX_SHADER = `#version 300 es
in vec3 aPosition;
in vec3 aNormal;
uniform mat4 uMVP;
out vec3 vNormal;
void main() {
    vNormal = aNormal;
    gl_Position = uMVP * vec4(aPosition, 1.0);
}`;

const FRAGMENT_SHADER = `#version 300 es
precision highp float;
in vec3 vNormal;
uniform vec3 uLight;
uniform vec3 uBase;
uniform float uFlat;
out vec4 outColour;
void main() {
    // Two-sided: the open fields show their inside, and a black backface reads
    // as a hole in the mesh rather than as the far wall of a cavity.
    float lit = abs(dot(normalize(vNormal), uLight)) * 0.8 + 0.2;
    outColour = vec4(uBase * mix(lit, 1.0, uFlat), 1.0);
}`;

/** Fail loudly, in the page, where the reader is. */
function fail(message) {
    const box = document.getElementById("lite-error");
    if (box) {
        box.textContent = message;
        box.hidden = false;
    }
    throw new Error(message);
}

// ─── the six lines of matrix maths ──────────────────────────────────────────

/** Column-major 4x4 product, `a` then `b`, as WebGL wants it. */
function multiply(a, b) {
    const out = new Float32Array(16);
    for (let c = 0; c < 4; c++) {
        for (let r = 0; r < 4; r++) {
            out[c * 4 + r] =
                a[r] * b[c * 4] +
                a[4 + r] * b[c * 4 + 1] +
                a[8 + r] * b[c * 4 + 2] +
                a[12 + r] * b[c * 4 + 3];
        }
    }
    return out;
}

/** Right-handed perspective, depth in `[-1, 1]`. */
function perspective(fovY, aspect, near, far) {
    const f = 1 / Math.tan(fovY / 2);
    const d = 1 / (near - far);
    return new Float32Array([
        f / aspect, 0, 0, 0,
        0, f, 0, 0,
        0, 0, (far + near) * d, -1,
        0, 0, 2 * far * near * d, 0,
    ]);
}

/** Right-handed look-at, with `up` fixed to `+y`. */
function lookAt(eye, centre) {
    const z = normalise([eye[0] - centre[0], eye[1] - centre[1], eye[2] - centre[2]]);
    // `up` and `z` are parallel only when the pitch is exactly +-90 degrees, and
    // the pitch clamp below stops just short of that.
    const x = normalise(cross([0, 1, 0], z));
    const y = cross(z, x);
    return new Float32Array([
        x[0], y[0], z[0], 0,
        x[1], y[1], z[1], 0,
        x[2], y[2], z[2], 0,
        -dot(x, eye), -dot(y, eye), -dot(z, eye), 1,
    ]);
}

const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
];
function normalise(v) {
    const l = Math.hypot(v[0], v[1], v[2]) || 1;
    return [v[0] / l, v[1] / l, v[2] / l];
}

// ─── the module ─────────────────────────────────────────────────────────────

const canvas = document.getElementById("isomesh-lite");
if (!canvas) {
    throw new Error("lite.js loaded on a page with no #isomesh-lite canvas");
}

const gl = canvas.getContext("webgl2", { antialias: true, depth: true });
if (!gl) {
    fail(
        "This demo needs WebGL2, which this browser has not enabled. The nine " +
        "playable demos need it too.",
    );
}

let wasm;
try {
    // No imports: the module is `no_std`-shaped, keeps no clock -- the browser
    // times the extraction -- and allocates through Rust's own allocator inside
    // its own linear memory.
    ({ instance: wasm } = await WebAssembly.instantiateStreaming(
        fetch("isomesh_web.wasm"),
        {},
    ));
} catch (error) {
    fail(`isomesh_web.wasm did not load: ${error}`);
}

const iso = wasm.exports;
const decoder = new TextDecoder();

/** One name out of the module's own table. `kind` 0 is fields, 1 extractors. */
function nameOf(kind, index) {
    const pointer = iso.iso_name(kind, index);
    const length = iso.iso_name_len(kind, index);
    if (pointer === 0 || length === 0) {
        fail(`isomesh_web has no name for ${kind}/${index}`);
    }
    return decoder.decode(new Uint8Array(iso.memory.buffer, pointer, length));
}

/** Fill a `<select>` from the module's name table, so there is one source. */
function fillSelect(select, kind, count) {
    for (let index = 0; index < count; index++) {
        const option = document.createElement("option");
        option.value = String(index);
        option.textContent = nameOf(kind, index);
        select.append(option);
    }
}

// ─── gl setup ───────────────────────────────────────────────────────────────

function compile(type, source) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        fail(`shader: ${gl.getShaderInfoLog(shader)}`);
    }
    return shader;
}

const program = gl.createProgram();
gl.attachShader(program, compile(gl.VERTEX_SHADER, VERTEX_SHADER));
gl.attachShader(program, compile(gl.FRAGMENT_SHADER, FRAGMENT_SHADER));
gl.linkProgram(program);
if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
    fail(`program: ${gl.getProgramInfoLog(program)}`);
}
gl.useProgram(program);

const uniform = {
    mvp: gl.getUniformLocation(program, "uMVP"),
    light: gl.getUniformLocation(program, "uLight"),
    base: gl.getUniformLocation(program, "uBase"),
    flat: gl.getUniformLocation(program, "uFlat"),
};
const attribute = {
    position: gl.getAttribLocation(program, "aPosition"),
    normal: gl.getAttribLocation(program, "aNormal"),
};

const buffers = {
    position: gl.createBuffer(),
    normal: gl.createBuffer(),
    index: gl.createBuffer(),
    wire: gl.createBuffer(),
};

gl.enable(gl.DEPTH_TEST);
// Backface culling stays off: `thin_plate` and `fbm_terrain` are open surfaces
// and culling them would draw half of each.
gl.disable(gl.CULL_FACE);
gl.clearColor(...CLEAR);

// ─── state ──────────────────────────────────────────────────────────────────

const controls = {
    field: document.getElementById("lite-field"),
    extractor: document.getElementById("lite-extractor"),
    samples: document.getElementById("lite-samples"),
    samplesLabel: document.getElementById("lite-samples-value"),
    wireframe: document.getElementById("lite-wireframe"),
};
const hud = {
    vertices: document.getElementById("lite-vertices"),
    triangles: document.getElementById("lite-triangles"),
    euler: document.getElementById("lite-euler"),
    nonManifold: document.getElementById("lite-non-manifold"),
    boundary: document.getElementById("lite-boundary"),
    degenerate: document.getElementById("lite-degenerate"),
    milliseconds: document.getElementById("lite-ms"),
};

fillSelect(controls.field, 0, iso.iso_field_count());
fillSelect(controls.extractor, 1, iso.iso_extractor_count());
controls.samples.value = String(SAMPLES_DEFAULT);

const camera = { yaw: 0.9, pitch: 0.45, distance: 6, centre: [0, 0, 0] };
let indexCount = 0;
let wireCount = 0;

/** Re-mesh, upload, and re-frame the camera if the domain changed. */
function remesh() {
    const samples = Number(controls.samples.value);
    controls.samplesLabel.textContent = `${samples}³`;

    const started = performance.now();
    const triangles = iso.iso_mesh(
        Number(controls.field.value),
        Number(controls.extractor.value),
        samples,
    );
    const elapsed = performance.now() - started;

    if (triangles === 0) {
        indexCount = 0;
        wireCount = 0;
        for (const span of Object.values(hud)) {
            span.textContent = "--";
        }
        draw();
        return;
    }

    // Views built here and nowhere else: `iso_mesh` may have grown the module's
    // memory, which detaches every view taken before this line.
    const vertices = iso.iso_vertex_count();
    const positions = new Float32Array(iso.memory.buffer, iso.iso_positions(), vertices * 3);
    const normals = new Float32Array(iso.memory.buffer, iso.iso_normals(), vertices * 3);
    const indices = new Uint32Array(iso.memory.buffer, iso.iso_indices(), iso.iso_index_count());
    const centre = new Float32Array(iso.memory.buffer, iso.iso_centre(), 3);

    gl.bindBuffer(gl.ARRAY_BUFFER, buffers.position);
    gl.bufferData(gl.ARRAY_BUFFER, positions, gl.DYNAMIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, buffers.normal);
    gl.bufferData(gl.ARRAY_BUFFER, normals, gl.DYNAMIC_DRAW);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, buffers.index);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, indices, gl.DYNAMIC_DRAW);
    indexCount = indices.length;

    // Three edges per triangle, undeduplicated: a shared edge drawn twice costs
    // one line and a de-duplicating pass over 40,000 triangles costs a frame.
    const wire = new Uint32Array(triangles * 6);
    for (let t = 0; t < triangles; t++) {
        const a = indices[t * 3];
        const b = indices[t * 3 + 1];
        const c = indices[t * 3 + 2];
        wire.set([a, b, b, c, c, a], t * 6);
    }
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, buffers.wire);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, wire, gl.DYNAMIC_DRAW);
    wireCount = wire.length;

    const extent = iso.iso_extent();
    camera.centre = [centre[0], centre[1], centre[2]];
    camera.distance = extent * 3;

    hud.vertices.textContent = vertices.toLocaleString();
    hud.triangles.textContent = triangles.toLocaleString();
    hud.euler.textContent = String(iso.iso_euler());
    hud.nonManifold.textContent = String(iso.iso_non_manifold_edges());
    hud.boundary.textContent = String(iso.iso_boundary_edges());
    hud.degenerate.textContent = String(iso.iso_degenerate_triangles());
    // The browser times this, which is why the module needs no clock and
    // therefore no `web-time` and no imports at all.
    hud.milliseconds.textContent = `${elapsed.toFixed(1)} ms`;

    draw();
}

function resize() {
    const ratio = Math.min(window.devicePixelRatio || 1, 2);
    const width = Math.round(canvas.clientWidth * ratio);
    const height = Math.round(canvas.clientHeight * ratio);
    if (canvas.width !== width || canvas.height !== height) {
        canvas.width = width;
        canvas.height = height;
    }
}

function draw() {
    resize();
    gl.viewport(0, 0, canvas.width, canvas.height);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    if (indexCount === 0) {
        return;
    }

    const cosPitch = Math.cos(camera.pitch);
    const eye = [
        camera.centre[0] + camera.distance * cosPitch * Math.cos(camera.yaw),
        camera.centre[1] + camera.distance * Math.sin(camera.pitch),
        camera.centre[2] + camera.distance * cosPitch * Math.sin(camera.yaw),
    ];
    const mvp = multiply(
        perspective(
            (45 * Math.PI) / 180,
            canvas.width / Math.max(canvas.height, 1),
            camera.distance * 0.01,
            camera.distance * 10,
        ),
        lookAt(eye, camera.centre),
    );

    gl.uniformMatrix4fv(uniform.mvp, false, mvp);
    gl.uniform3fv(uniform.light, normalise(LIGHT));

    gl.bindBuffer(gl.ARRAY_BUFFER, buffers.position);
    gl.enableVertexAttribArray(attribute.position);
    gl.vertexAttribPointer(attribute.position, 3, gl.FLOAT, false, 0, 0);
    gl.bindBuffer(gl.ARRAY_BUFFER, buffers.normal);
    gl.enableVertexAttribArray(attribute.normal);
    gl.vertexAttribPointer(attribute.normal, 3, gl.FLOAT, false, 0, 0);

    gl.uniform3fv(uniform.base, SURFACE);
    gl.uniform1f(uniform.flat, 0);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, buffers.index);
    // `UNSIGNED_INT` is the reason this is WebGL2: a 49^3 Marching Tetrahedra
    // mesh is well past what a `Uint16Array` index can address.
    gl.drawElements(gl.TRIANGLES, indexCount, gl.UNSIGNED_INT, 0);

    if (controls.wireframe.checked && wireCount > 0) {
        gl.uniform3fv(uniform.base, WIRE);
        gl.uniform1f(uniform.flat, 1);
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, buffers.wire);
        gl.drawElements(gl.LINES, wireCount, gl.UNSIGNED_INT, 0);
    }
}

// ─── input ──────────────────────────────────────────────────────────────────

let dragging = false;
canvas.addEventListener("pointerdown", (event) => {
    dragging = true;
    canvas.setPointerCapture(event.pointerId);
});
canvas.addEventListener("pointerup", (event) => {
    dragging = false;
    canvas.releasePointerCapture(event.pointerId);
});
canvas.addEventListener("pointermove", (event) => {
    if (!dragging) {
        return;
    }
    camera.yaw -= event.movementX * 0.008;
    // Stopping just short of the pole keeps `lookAt`'s `up` from being parallel
    // to the view direction, where the cross product collapses.
    const limit = Math.PI / 2 - 0.02;
    camera.pitch = Math.min(limit, Math.max(-limit, camera.pitch - event.movementY * 0.008));
    draw();
});
canvas.addEventListener(
    "wheel",
    (event) => {
        event.preventDefault();
        camera.distance = Math.min(
            iso.iso_extent() * 20,
            Math.max(iso.iso_extent() * 0.2, camera.distance * Math.exp(event.deltaY * 0.001)),
        );
        draw();
    },
    { passive: false },
);

controls.field.addEventListener("change", remesh);
controls.extractor.addEventListener("change", remesh);
controls.samples.addEventListener("input", remesh);
controls.wireframe.addEventListener("change", draw);
window.addEventListener("resize", draw);

remesh();
