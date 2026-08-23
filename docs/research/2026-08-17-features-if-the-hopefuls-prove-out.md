# Features, if the hopefuls prove out

**Date:** 2026-08-17
**Companion to:** `2026-08-17-mechanics-from-the-field.md`
**What this is:** the nine candidates translated into game features — verbs, UI, progression, failure
states — assuming each one's pre-registered predictions come back positive. Organised by *feature*
rather than by candidate, because most features need two proofs and one needs four.

---

# Part 0 — The thesis: ten designer-placed boxes become ten measured numbers

Every candidate in the dossier does the same thing, and it is worth naming once because it explains why
they matter to *this* game and not to games in general.

| What games place by hand today | What replaces it | Needs |
|---|---|---|
| Reverb volume | the connected air component | §1.3 |
| Traversal gate ("you need the small form") | clearance in metres at the throat | §1.1 |
| Tactical point / cover marker | a local minimum of `ρ` on the skeleton | §2.4 |
| Breakable prop, scripted collapse | an admissible force state exists, or does not | §2.1 |
| Fire/blight spread radius | geodesic distance on the actual surface | §2.2 |
| Decal projector, UV atlas | the logarithmic map | §2.2 |
| "Aged rock" texture variant | reaction–diffusion state keyed to the grid edge | §2.3 |
| Pre-fractured asset | the shape's own fracture modes | §3 |
| Hand-authored cave layout | dissolution along the flow the player routed | §1.2 |
| Scattered ambient formations | speleothems where water actually drips | §2.5 |

**A hand-placed box is exactly the thing destruction invalidates.** That is why a destructible game
either ships shallow systems or ships a lot of designer labour that the player is then quietly
prevented from breaking. Every row above survives the player digging through it, because it is
recomputed rather than remembered. That is the whole argument, and it is why these are features rather
than tech demos.

---

# Part 1 — Features from a single proof

## 1.1 CLEARANCE — a number on your reticle, and gating with no gates

**Proof needed:** §1.1's identity holds discretely (the O(h) oracle check) and λ-membership is
incremental under digging.

**The verb.** Look at a passage; the HUD reads **0.62 m**. That is `2r` at the tightest point of what
you are looking at — the diameter of the largest ball that fits. Your character card says 0.55 m. It
fits. The thing chasing you says 0.90 m. It does not.

**What that makes possible, in order of how much it changes the game:**

**(a) Metroidvania gating with zero authored gates — and the player can build the gates.** A traditional
progression gate is a door a designer placed that needs the double jump. Here, *every passage narrower
than your λ is a gate, everywhere, automatically, including passages the player dug themselves.* Acquire
a smaller form and **the entire existing map re-opens at once** — every 0.5 m crack you have walked past
for ten hours becomes a route, and the player knows exactly which ones because they have been reading
the number all along. That inverts the usual reveal: instead of "where is the door this key opens", it
is "I already know six places I could not fit, and now I can."

And it runs backwards: **you make gates for other things.** Narrow a passage to 0.8 m and the brute is
permanently locked out of your base while your 0.5 m allies come and go. Defensive architecture becomes
a measurement problem instead of a hitpoint problem.

**(b) Combat you win with geometry.** Escape is not a stamina bar, it is a fit. Pursuit is not damage,
it is widening. A fight where the brute is chewing the tunnel wider while you are deeper in and
counting down its progress is a tension curve nothing currently produces, because it is *legible* —
the player can see the number climbing.

**(c) Logistics with clearance ratings.** Which cargo fits which route. A hauler that is 1.4 m needs a
1.4 m corridor for its entire path, so a base's internal roads have a spec, and one careless cave-in
downgrades a main line. This is the mechanic that makes underground infrastructure worth maintaining.

**(d) Bestiary entries with a real number in them.** "Gnasher — 0.90 m." Players will learn creature
widths the way they currently learn damage types, and the knowledge transfers to every situation in the
game rather than to a lookup table.

**Risk, stated plainly.** This feature is **invisible without UI**. The geometry does not look different;
only the rules changed. If the clearance readout is not on screen and legible at a glance, players will
experience it as arbitrary. The HUD element is not polish here — it *is* the feature.

## 1.2 HYDROLOGY — water routing as a construction verb, and a map that changes while you are away

**Proof needed:** §1.2's headless graph simulator reproduces wormhole competition and the bimodal
aperture distribution.

**The verb.** You do not dig the cave. You **decide where the water goes**, and the cave digs itself
along that route over the next several minutes of play.

**What that makes possible:**

**(a) A long-horizon goal that is not combat and not a crafting tree.** "Get the water running through
the east massif so it opens a route to the vault" is a multi-session objective involving surveying,
damming, channel-cutting, and waiting — and the payoff is a passage that did not exist and that nobody
placed. It is the first goal in a voxel game that the *world* completes rather than the player.

**(b) Competition as tension.** The literature's result is that one wormhole wins and the rest go dry.
So **your channel competes with natural fractures, and you can lose.** You cut a feeder, it starts
widening, and then a fracture three hundred metres east captures the flow and your passage stalls
forever as a dead-end crawl. That failure is far more interesting than "the dig took too long", and it
teaches hydrology: you learn to check what else is drawing from the same head.

**(c) Sabotage and denial in PvP.** Dam the enemy's supply and their access route stops growing. Redirect
it into their base and it dissolves *toward* them. Attacking someone's infrastructure by changing where
water goes is a category of aggression no game has.

**(d) A server whose map is genuinely different at month six.** Passages nobody dug. Chambers that opened
while everyone was offline. And because it composes with air connectivity, **the world can breach itself
while you are away** — dissolution connects two regions, the breakthrough event fires, and something got
in through a hole no player made.

**(e) Restoration as a goal.** Flow can be given back. A dry passage can be re-wetted. That turns
hydrology into something you can repair, which is rarer and better than something you can only exploit.

**Risk.** Timescale. Ambient dissolution must be tuned to essentially zero, or the world reads as
weather rather than consequence. The feature only works if it is *provoked* — visible widening within
10–30 seconds of a deliberate action, breakthrough in 5–20 minutes.

## 1.3 THE ROOM YOU DUG — hearing space, and echolocation as navigation

**Proof needed:** §1.3's Sabine RT60 per air component in under 0.1 ms; a 2D re-bake under 30 ms.

**The verb.** Throw a stone into the dark and **listen to how big the space is.**

**What that makes possible:**

**(a) Navigation without light.** A genuine blind-traversal mode, or a light-scarce game where sound is
the primary survey instrument. Currently impossible because reverb is a designer volume — you cannot
echolocate a space nobody has authored a reverb for, and in a dug world that is every space.

**(b) Spoiler-free spatial hints.** A long tail means a big chamber ahead. Players learn to hear the
difference between a dead-end and a cavern and route accordingly, which makes exploration decisions
informed without a map marker.

**(c) Confirmation feedback for sealing, with no UI at all.** You wall off a chamber and **the room goes
dead.** That is the airtightness mechanic reporting success through the audio, on the frame it happens.
It is the cheapest possible interface for the most abstract mechanic in the game.

**(d) Tap-to-survey.** Strike a wall; a void behind it changes the response before you break through.
Ore and cavity detection that is *physical* rather than a highlighted overlay — and therefore a skill
rather than an unlock. (This one wants the modal work in §3 as well, so it is a Part 2 feature
properly; the air-side version — shout and listen for a return — works on §1.3 alone.)

**(e) Multiplayer shared signal.** The collapse sounds different to the person in the next chamber than
to the person standing in it, correctly, without anyone authoring either.

**Why this ships first.** It is the cheapest of the nine, it has a public reference implementation, it
monetises infrastructure that is already built and measured, and it is the one a player perceives
immediately without being told to look for it.

## 1.4 THE SURVEY SCANNER — a field that says what will fall

**Proof needed:** §1.4's top-1% capture rate ≥ 80% against brute-force ablation.

**The verb.** Raise the scanner and the rock face **paints**: cool where removal is harmless, hot where
removal moves the structure toward collapse.

**What that makes possible:**

**(a) Mining as a read-then-act loop.** The interesting decision in a mining game is currently "where is
the ore." This adds "and can I take it from here." Two competing gradients over the same wall is a real
choice.

**(b) A skill ceiling with a legible ladder.** Early: you use the scanner constantly. Later: you read the
thrust veins directly (§2.1's diegetic visualisation) and only scan when unsure. Expert: you know from
the shape of the room. Nothing is unlocked; the player just stops needing the crutch — which is the
best kind of progression because it is real.

**(c) Precision sabotage in PvP.** Find the single cell whose removal drops the enemy's ceiling. That is
a skilful, physics-true action with none of the "shoot the wall until its HP is zero" texture. And it is
*defendable*, because the defender can run the same scanner and shore up what it shows.

**(d) A tutorial that teaches itself.** The field is the explanation. A player who watches it change as
they carve learns load paths without a tooltip.

**Risk.** This is a UI, and a noisy or ugly one is worse than nothing. It also must not be so
authoritative that it removes the decision — the interesting version says *how much*, not *yes/no*.

---

# Part 2 — Features that need two or more proofs

## 2.1 LOAD-BEARING ARCHITECTURE — the ceiling falls because it was holding something up

**Proofs needed:** §2.1's admissibility solve at ≤ 250 ms for 1000 macro-blocks, **and** the arch golden
values, **and** the legibility gate (players above ~70% correct prediction), which in practice means
§1.4 and the thrust visualisation too.

**The verbs.** Shore. Pillar. Undercut. Buttress. Bench.

**What that makes possible:**

**(a) Room-and-pillar mining as a strategy the game actually rewards.** Leave pillars and the roof
holds; take them and it does not. This is a real technique with a real name that no game has ever needed
the player to know, because no game has computed the thing it addresses.

**(b) Collapse as a scenario generator rather than an authored setpiece.** A cave-in traps you — and now
the objective is to dig out, from the inside, with the rock's remaining structure as the constraint. The
game did not place that. Your dig did.

**(c) Material identity from two parameters.** Friction coefficient and density are the *only* physical
inputs, so "wet clay" and "granite" differ in a way a player can hold in their head and predict. That is
unusually rare — compare a stiffness-based system whose behaviour depends on solver iteration count and
is therefore unlearnable.

**(d) Buildings that stand for a reason.** The same solver grades player construction. A vault holds
because it is a vault. And the failure mode is instructive rather than arbitrary: the solver names the
hinge lines, so the collapse *shows you where you were wrong*.

**(e) The funicular build tool** (§2.1's Candidate E, the cheapest thing in the whole dossier at
sub-millisecond). Drag a span and it **snaps toward the shape that stands**. The player learns arch form
by being gently corrected a hundred times, and the entire family of standing arches is one slider — the
force-density scale. That is a build tool that teaches structural intuition without a single word of
text.

**Risk, and it is the real one.** Correct and unpredictable is worse than approximate and legible. The
harness has a literal player-model test for this — show twenty people a pillar, ask whether cutting here
drops it — and **below ~70% the feature is unshippable in that form.** The fix is not to soften the
physics; it is §1.4's field and the thrust veins, which is why these three are one feature and not
three.

## 2.2 SURFACES THAT REMEMBER — age, spread, and forensics

**Proofs needed:** §2.3's grid-edge keying survives 1,000 edits with <3% spurious resets, **and** §2.2's
incremental factorisation (or the CPM route) for the geodesic half.

**What that makes possible:**

**(a) Rock that shows its age, with no timestamps and no UI.** A wall cut yesterday carries a grown
pattern; a wall cut a minute ago is bare; where they meet there is a **visible merge seam** as the old
growth advances into the new cut. A chamber's excavation history is written on it.

**(b) Forensics and scouting as a genuine read.** Walk into someone else's workings and *see* which parts
are old and which are from this morning. In PvP that is intelligence you gather by looking rather than
by a detection ability. In a single-player game it is environmental storytelling that generates itself.

**(c) Spread mechanics with the correct metric.** Fire, blight, frost, corruption that travel *along the
surface* at a fixed speed — so they climb walls, cross ceilings, and pour down into your tunnel. And the
counter-verb is real: **cut a firebreak**, and you can see in the moment of cutting whether it is wide
enough, because the front slows and detours. Sever the last land bridge and it stops dead while sitting
two metres away in straight-line space. Nothing ships this; Minecraft-style fire uses volumetric
adjacency and cannot tell "around a ledge" from "through a thin wall."

**(d) Territory that grows correctly.** Moss, fungus, a faction's claim — spreading by geodesic distance
means it wraps corners the way a real growth would, and its extent is a *measurable area* rather than a
radius.

**(e) Decals with no UV atlas** (from the log map). Graffiti, blast scorches, blood, paint — applied to
carved geometry that had no texture coordinates a second ago.

**(f) Carving history as pattern morphology.** The dossier's neatest detail: the edit log gives a
per-vertex time-since-exposure for free, which is exactly the spatially varying diffusion coefficient
that widens or narrows the pattern. **How you carved determines how the rock ages, at zero storage.**

## 2.3 THE THROAT — tactics that expire when you breach

**Proofs needed:** §2.4's throat diameter within ±1 voxel, **and** the noise filter removing ≥90% of
false throats.

**What that makes possible:**

**(a) Siege maths that both sides can compute.** A throat is `2ρ` metres wide, so "how many abreast" is
arithmetic. A 1-wide throat gets one defender who cannot be flanked. A 3-wide throat gets a line. Both
attacker and defender are looking at the same number.

**(b) Attacking a position by changing its geometry.** Blow a second hole in the flank and the old
throat's importance collapses — **the defensive position becomes invalid automatically**, and the AI
relocates because the derived graph changed, not because a scripted volume was invalidated. That is the
first tactical AI that survives a destructible world without designer maintenance.

**(c) Defence as excavation.** You do not place cover; you *cut* it. Narrowing an approach is a
construction action with a measurable tactical result.

## 2.4 THE CARVED BELL — geometry you can hear

**Proofs needed:** §3's 48 modes under 8 ms, **and** the unresolved perceptual question of how few modes
sound like a bell, **and** a single dig actually perturbing λ₁ above the just-noticeable difference.

**What that makes possible, if the perceptual gate passes:**

**(a) Prospecting by ear.** Tap the wall; a void behind changes the response before you break through.
Combined with §1.3, that is a two-instrument survey — strike for what is *behind* the rock, shout for
what is *beyond* the opening.

**(b) Carving as instrument-making.** Hollow a slab and its pitch drops. Thin a tongue and a partial
appears. Some players will make bells, and they will be *correct* bells. This is the single most likely
thing in the document to end up in a video.

**(c) Diagnostics.** A structure under load sounds different from one that is not. Creaking that means
something.

**Honest position.** The cost analysis lands at ~40 modes, and whether 40 modes reads as "a bell" or as
"a filtered noise burst" is unresolved in the literature I could reach. **Run the listening test before
the eigensolver.**

## 2.5 A CAVE THAT AGES — speleothems, and conservation as a mechanic

**Proofs needed:** §2.5's one-parameter profile within 5% of the integrated growth law, plus §1.2's
hydrology to route the water.

**What that makes possible:**

**(a) Time made visible.** A chamber you hollowed out early in a save becomes visibly *old* by the end
of it. That is a form of long-session payoff that costs nothing to run and cannot be faked with a
texture swap, because the formations are where the water actually goes.

**(b) A renewable resource with a real conservation rule.** Harvest formations, and they regrow **only
where water still drips.** Over-harvest a chamber and it stays bare — unless you restore the hydrology
that fed it. That is an ecology mechanic derived from physics rather than a respawn timer.

**(c) Dating a space.** Formation size tells you roughly how long this chamber has been open, which
pairs with §2.2's surface ageing into an actual forensic vocabulary.

---

# Part 3 — What each is worth, honestly

| Feature | Player-perceptible? | Cost to prove | Cost to ship | Risk |
|---|---|---|---|---|
| **The room you dug** (§1.3) | **Immediately, without being told** | lowest | low | almost none |
| **Clearance** (§1.1) | Only via UI | low | low + UI design | invisible if the HUD is wrong |
| **Surfaces that remember** (§2.2) | Yes, gradually | medium | medium | needs the tufted Laplacian or the chemistry diverges |
| **Hydrology** (§1.2) | Yes, if provoked | low | **high** — needs a whole water system | reads as weather if mistuned |
| **Survey scanner** (§1.4) | Yes | medium | low | a bad UI is worse than nothing |
| **Load-bearing architecture** (§2.1) | Yes, dramatically | high | high | correct-but-illegible is a real failure state |
| **The throat** (§2.4) | Indirectly, via AI behaviour | low | medium | false throats from destruction noise |
| **The carved bell** (§2.5/§3) | Yes, memorably | high | high | 40 modes may not sound like anything |
| **A cave that ages** (§2.5) | Slowly | low | low | needs hydrology first |

**Two features are nearly free and immediately felt: the acoustics and the clearance readout.** Two are
the headline: load-bearing architecture and hydrology. The rest are depth.

---

# Part 4 — The games these compose into

None of these is a game. Four combinations are.

**The cave ecology game.** Hydrology + speleothems + acoustics + connectivity. You are a steward of a
living cave system: you route water, the system responds over sessions, formations grow where you let
them, chambers open and go dry, and you hear all of it. The core loop is *tending* rather than
extracting, which is a genre voxel games have never had because their worlds do not do anything on
their own.

**The λ-Metroidvania.** Clearance + throat + fill. Traversal gating with no authored gates, in a world
where the player creates and destroys their own gates and where acquiring a different body form
re-opens everything at once. The map is the ability system.

**The mining engineering sim.** Load-bearing architecture + survey scanner + tap-to-survey + clearance.
*Deep Rock Galactic* with real geotechnics: you shore, you bench, you leave pillars, you read the rock,
and the mountain kills you when you are wrong. The failure states are educational and the skill ceiling
is genuine knowledge.

**Asymmetric siege.** Throat + admissibility + clearance + acoustics. The attacker's job is to change
the defender's geometry — widen an approach, drop a ceiling, breach a flank — and the defender's job is
to narrow, shore and seal. Both sides are doing structural engineering under time pressure, and both
can compute the same numbers.

---

# Part 5 — The one-paragraph version

**If all nine prove out, the game's rules stop living in designer-placed boxes and start living in the
world's own measured geometry** — which is the only version that survives a player who can dig through
anything. The two that pay immediately are hearing the size of a space on the frame you break into it,
and a clearance number on your reticle that turns creature size into a traversal rule and makes every
narrow crack in the world a gate you already know about. The two that are the headline are a world that
digs its own caves along the water you routed, and rock that falls because it was holding something up
and will tell you, before you cut, how much closer this cut moves it. And the one that will end up in a
video is a player carving a slab until it rings on the note they wanted.
