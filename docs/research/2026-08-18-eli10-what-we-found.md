# What we found, explained like you're ten

This is every big discovery the `isomesh` project has made so far, written so a ten-year-old can follow it — and always answering the question that actually matters: **what does this mean when you're playing the game?**

The game is one where the whole world is made of stuff you can dig through. Not blocks — actual smooth rock, with caves and tunnels and overhangs. You can dig anywhere, and the ground has to redraw itself instantly.

---

## First, the one idea everything else sits on

The computer does not store the shape of the world. That would be enormous.

Instead it stores a **number at every point in space**, and the number answers one question: *how far are you from the nearest rock, and are you inside it or outside it?* Negative means you're inside the rock. Positive means you're in the air. Zero means you're standing exactly on the surface.

So the world is a giant invisible field of numbers, and **the surface of the world is just "everywhere the number is exactly zero."**

The job of `isomesh` is to look at that cloud of numbers and draw the triangles that sit where the zeros are. That's it. Numbers in, triangles out. Everything below is what we learned trying to do that fast, correctly, and in ways that make the game more fun.

---

## Discovery 1 — The slow part wasn't the hard part. It was two silly things.

We had two families of algorithms. One of them, the "dual" family, was famous for being slow. Everyone knows this. It's in books.

We measured it, and the slowness came from two things that have nothing to do with the algorithm at all.

**Silly thing one.** Deep inside the code there was a loop that had to keep asking "am I working on the X direction, the Y direction, or the Z direction?" — *every single time around*, millions of times. We changed it so the computer is told which direction ahead of time and makes three separate copies of the loop, one per direction. Nothing else changed. **It got about three times faster, and produced byte-for-byte the exact same triangles.**

**Silly thing two.** Computers keep a small fast pocket of memory near the processor. Our rows of world data were spaced *exactly* the wrong distance apart — every row landed in the same slot in that pocket, so each one kicked the last one out. It's like a class where all thirty kids are assigned the same locker. The fix was to make each row one number longer, so they stop colliding. **That gave us another big chunk of speed, and cost less than one percent of memory.**

Put together: the "slow" algorithm went from **694 milliseconds to 163 milliseconds** on a big chunk of world. It went from five and a half times slower than its rival to only about a quarter slower — and below a certain size, it's now the *faster* one.

**What it means when you're playing:** the world can be bigger, the chunks can rebuild faster after you dig, and the smoother-looking mesher stopped being a luxury you couldn't afford. "Dual contouring is slow" turned out to be a statement about one piece of code on one afternoon, not about the algorithm.

---

## Discovery 2 — Our ruler was broken, and it made two true things look false.

We wanted to check whether our triangles pointed the right way along a sharp edge. So we built a test shape — a wedge, like a doorstop — where we could *calculate* the correct answer by hand and compare.

Twice in a row, the test said our mesher was wrong.

It wasn't. **The thing we were comparing against was broken.** The code that computed the "true" direction had a spot where two nearly-equal numbers got subtracted, leaving basically nothing but rounding dust — and then it confidently pointed in whatever random direction that dust happened to face. For about half the points near a flat face, our "ground truth" was noise.

When we fixed the ruler, the count of bad vertices dropped from **6,959 to 472**, and both hypotheses that had looked dead turned out to be **exactly right**.

There's an earlier version of the same lesson. We built a test wedge lined up neatly with the grid, and it agreed with our prediction to four decimal places. That felt great. It was a trap — the shape was *too* well behaved. When we rotated it by 17 degrees, it fell apart immediately.

**The rule we wrote down:** if a measurement seems impossible, suspect the measuring stick before you suspect the world. And if a test agrees *perfectly*, be suspicious — perfect agreement usually means the test isn't testing anything.

**What it means when you're playing:** nothing directly — and that's the point. This is the discovery that stopped us shipping a "fix" for a bug that didn't exist, which would have made the actual sharp edges in the game worse.

---

## Discovery 3 — Digging is cheap. Writing down the answer is expensive.

When you dig a small hole, only a tiny bit of the world actually changes. We measured it: an edit that touches **0.038%** of the world changes about 792 little cells — and that number stays *exactly the same* whether the world is small or 64 times bigger. That's great news. The work is proportional to the digging, not to the world.

But then we looked at the list of triangle corners the mesher hands to the graphics card, and **56 to 77 percent of the entries had moved.**

Why? Because the corners are numbered 1, 2, 3, 4... in the order they're found. If one cell suddenly produces one extra corner, **every single number after it shifts by one.** Nothing about that geometry changed. Only its place in the queue did.

We worked out the ceiling: if we named each corner after *where it lives in the world* instead of *when it was found*, the churn would drop from **15,706 to 346** — about 45 times better, and it would stay flat no matter how big the world gets.

Then we hit a wall, and the wall is interesting. The only version of that idea that actually works requires the mesher to *remember things between calls* — and we have a test, written on purpose, that exists specifically to catch a mesher whose answer depends on what it did last time. So the fix would turn a safety test's failure condition into its intended behaviour. We stopped rather than paper over it.

**What it means when you're playing:** right now, digging a small hole makes the game send far more data to the graphics card than it needs to. Fixing it means bigger destructible worlds at the same frame rate. We know exactly how much is on the table, and we know exactly what's blocking it.

---

## Discovery 4 — Digging breaks caves in two more often than you'd think

We asked: when you fill in a bit of a cave, how often does that split the cave into two separate caves?

**About one time in six.** And splits outnumber "you sealed off a tiny pocket entirely" by 27 to 5.

That matters because a lot of interesting game logic depends on knowing which spaces connect to which. Can the water flow from here to there? Can the monsters path from their nest to your base? Did you just seal yourself in?

Checking this the dumb way — recompute everything — would be far too slow to do every time you dig. We built the smart version instead, and it does about **436 times less work** than recomputing.

But then we built the *nastiest possible* test: an edit that slices a cave exactly down the middle into two equal halves. That's the worst case, because you have to explore both halves before you can be sure they're really separate. It cost about **1.1 times a full recompute** — slightly worse than just starting over.

The fix wasn't a cleverer algorithm. It was **chopping the world into chunks**, so any search is trapped inside one chunk and can't run away across the whole world.

And one more honest note: that 1.1× number is the *adversarial* case, deliberately built to be awful. In a realistic scene — severing a passage between two chambers — the same measurement came out at **0.028×**, which is 35 times cheaper. We nearly wrote down the nasty number as if it were the normal one. We now have a rule about that: run the adversarial fixture *and* a realistic one, and never quote one number for both.

**What it means when you're playing:** the game can know, on the exact frame it happens, that you just cut a tunnel in half — or that you just broke through into a chamber nobody has entered. Instantly, every time, without a hitch in the frame rate.

---

## Discovery 5 — The echo is free

Once the game knows which pocket of air you're standing in, it also knows two things it was already tracking anyway: how big that pocket is, and how much wall surface it has.

There's a very old, very simple formula in acoustics — Sabine's formula — that turns exactly those two numbers into **how long an echo lasts in a room**.

So we measured how expensive it is to compute the reverb on the frame you smash through into a new cavern.

**0.0015 microseconds.** That is two numbers read out of memory and one division. The budget we'd allowed was 100 microseconds. We came in about sixty thousand times under it.

And the answer it gives is sensible: two chambers merging produced a **0.70 second** reverb tail, which is a believable big-cave sound.

**What it means when you're playing:** the moment your pickaxe punches through the last bit of rock into a huge chamber, the sound of your own footsteps changes — instantly, correctly, calculated from the actual shape of the actual hole you actually just made. Not a trigger volume an artist placed. Not a preset. The real room. Dig a side passage and the echo shortens. Seal it up and it changes back.

---

## Discovery 6 — You cannot hear someone digging inside a pillar, and we were wrong about why

Here was a lovely idea: if the world can be "played" like an instrument, a stone pillar has a note it rings at. Dig material out of it and the note should change. So you could *hear* structural damage.

We tested it. Digging inside a pillar changed its note by **0.11%**, and a human ear can only notice about **0.6%**. So: inaudible. Fine, that half we expected.

The half we *didn't* expect: we were sure that digging near a thin bit — the narrow web connecting two thick parts — would be dramatic, because thin bits look fragile. It changed the note by **0.032%**. Even less. The prediction that died was the one we were confident about.

The reason is genuinely surprising. What decides whether removing a chunk changes the note is not *whether the chunk is in a skinny place*. It's **how much of the object's bending energy is stored right there**. The thin web isn't doing the bending — it's just going along for the ride while the thick pillar sways. Taking material out of a passenger doesn't change the ride.

And we checked our instrument wasn't just numb: carving out a full 20% cavity moved the note **63.9%**. It can absolutely detect a change. There just wasn't one.

**What it means when you're playing:** "hollow it out and hear it groan" doesn't work, and we found that out in a day instead of a month. What survives is more interesting: the note *is* a real measure of structural health, and the places where bending energy concentrates are exactly the hinges where a structure will actually fail. That's a creak that means something.

---

## Discovery 7 — We can tell you when the arch falls, and we checked it against real books

If you're going to let players build stone arches and bridges, you need to know when one collapses. Not "it looks wobbly" — actually collapses.

There's real engineering literature on this, going back centuries. We built the solver and checked it against two published numbers.

- The thinnest an arch can be before it falls: books say **0.1075**. We got **0.10734**.
- How far you can tilt an arch before it collapses: books say **15.84 degrees**. We got **15.850 degrees**.

We also learned something in the *middle* of building it: our first solver used a fashionable fast method that gets close to the answer quickly and then crawls. Crawling wasn't good enough here — we needed to land precisely on the line between "stands" and "falls." We swapped it for a simpler, older method that actually converges. Good thing we noticed before trusting the numbers.

**What it means when you're playing:** you dig out the base of a stone bridge. The game *knows* — using the same maths a real engineer would use — whether it still stands. Take one more block and it comes down. Not scripted. Not a hit-point bar. Actual statics.

---

## Discovery 8 — The shape that passed every test and was still broken

We have a whole battery of tests for whether a mesh is a proper solid object: is it closed, do the triangles agree on which way is out, does its Euler characteristic come out right, does every edge have exactly two faces.

We built a shape called a **bowtie** — two pyramids touching at a single sharp point, tip to tip. Forty-eight numbers, that's all.

It passed **every single test**. And it is not a valid solid, because that shared tip is a place where the surface pinches to nothing — you can't say what's "inside" there.

The tell was hiding in plain sight: its Euler characteristic came out as **3**. That number must be *even* for any closed solid. It was odd. The maths was screaming and nothing was listening, because no test checked that particular thing.

The fix was one line: pass a number that was already being computed, and require it to be zero.

**What it means when you're playing:** physics engines take meshes and decide what's solid. Hand one a bowtie and you get a player who falls through the floor, or a bullet that passes through a wall, in a way nobody can reproduce. This is the bug class that eats a week. Now it can't get out of the building.

---

## Discovery 9 — The documentation said the GPU was slower. Its own data said 37 times faster.

Our README answered the question "is doing this on the graphics card faster?" with **no**.

The project's own committed measurement file, sitting right there in the same repository, said the graphics card was **37.6 times faster** at a decent grid size.

The README was right *once*, before two improvements landed, and nobody re-read it. Worse, that wrong sentence got copied out to the public documentation page the same day.

There's a related one. A result got fixed, the fix changed what a whole comparison meant, and the explanation was written **only in a commit message**. Six weeks later nobody could find it. The rule now: **a result that exists only in a commit message is not retrievable.** It goes in the findings file or it didn't happen.

**What it means when you're playing:** nothing directly, but it's why this project keeps a findings file at all. Documentation rots quietly, and the version that rots is always the one that used to be true.

---

## Discovery 10 — We tested on a real bonsai tree

Most graphics papers test on maths shapes — spheres, donuts, mathematical squiggles. Nice and clean. Suspiciously clean.

So we went and got real scientific CT scan data: a scanned **bonsai tree** and a scanned **jet of fuel**. Real, messy, noisy, quantised measurements.

Two things came out of it.

First, a published claim about triangle quality landed **inside** the range the literature reports — which told us that range is a property of *the kind of data*, not of anyone's clever implementation. An earlier test had seemed to contradict it; the earlier test had just used the wrong kind of shape.

Second, on a surface with **over a million triangles**, plain old Marching Cubes produced **zero** bad edges. And the fancier "guaranteed-manifold" algorithm proved it earns its keep on real data specifically: it took bad edges from **1,776 down to 85**, a 95% cut — which the clean maths shapes had never shown, because they never had the problem.

We also found that a third algorithm flatly **refused** to mesh the bonsai. Chasing that down turned up something worth knowing: the scan stores values as whole numbers, and at a sharp peak the neighbours on both sides are *equal*, so the usual way of computing a slope gives exactly zero — and you can't point "outward" when there's no outward. **16,284 sample points sat exactly on the surface.** The workaround is charmingly simple: cut the surface at 127.5 instead of 127, so no sample can land exactly on it.

**What it means when you're playing:** the game can eat real-world scan data, not just procedurally generated shapes. And we know which algorithm to reach for when the data is ugly, because we tried it on ugly data.

---

## Discovery 11 — Two things we proved were impossible, and were glad to

Not every result is "we made it work." Some are "we saved ourselves months."

**The medial axis detector.** There's an elegant idea: the "spine" running down the middle of a shape is where a certain quantity goes exactly to zero, so just look for the zeros. We found **4,225** of them, exactly where the spine should be. Beautiful.

Then we shifted the grid by half a step. **All 4,225 vanished.** We shifted by a full step and got exactly **3,136** back — which was the precise number our own maths had predicted beforehand for that alignment. The detector wasn't finding the shape's spine. It was finding **where the grid happened to be sitting relative to the origin**. It was measuring itself.

**Better vertex placement.** A proposal claimed that placing triangle corners more cleverly would improve accuracy by more than 20%. We killed it with arithmetic, before writing a line of code: even if you place *every single corner perfectly* on the true surface, a flat triangle still sags away from a curved surface in its middle. That sag is a floor you cannot get under. The proposed 20% was **above the ceiling of any placement rule whatsoever**.

And a sharper version of the same fact fell out: one family of algorithms is limited by *where its corners are* (so better placement helps), and the other is limited by *how it cuts the surface into triangles* (so better placement can't help at all).

**What it means when you're playing:** nothing — which is exactly the return on investment. Both of these would have been weeks of work producing something that could not have worked, and both died in under a day.

---

## The habits, and why they're the real result

Reading all of the above back, the pattern is that **most of the value came from being wrong efficiently.**

So here are the rules the project earned the hard way, each one paid for by a specific mistake:

- **Write down what you expect *before* you measure it**, including the exact result that would prove you wrong. Otherwise you'll find a way to feel right afterwards.
- **Never edit a prediction after seeing the answer.** That's not a rule about honesty, it's a rule about usefulness — an edited prediction can't teach you anything.
- **Count things instead of timing them where you can.** Counts are exact and mean the same thing on every machine. A committed timing file turned out to be **1.45× stale** without anyone noticing — and separately, *adding one unrelated function* moved a measurement from **152.5 ms to 130.8 ms**, purely by shifting the code around in memory. Neither of those is the algorithm changing.
- **Check the ruler as hard as the thing you're measuring.**
- **A paywalled paper is not a missing paper.** Check the author's own page — we found one sitting there under the exact filename we already knew.
- **A rule written down once has not been applied once.** We wrote a rule, then broke it on the very next experiment. Now the rules are tests.
- **Never delete a result that turned out to be wrong.** It records which sources to distrust, which is worth more than the fact was.

That last one is why the findings file only grows. It's not a log of successes. It's a map of every place the ground turned out to be soft.

---

*Companion documents: the corpus audit, the novelty table, the backlog and harness audit, and the memo on the four blocked decisions. All dated 2026-08-18.*
