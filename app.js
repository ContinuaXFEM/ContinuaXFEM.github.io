// ContinuaX website: the live earthquake hero, the replayed agent session and scroll reveals.
"use strict";

const reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

// ---------------------------------------------------------------- reveal on scroll
(() => {
  const els = document.querySelectorAll(".reveal");
  if (!("IntersectionObserver" in window)) { els.forEach(e => e.classList.add("in")); return; }
  const io = new IntersectionObserver(es => es.forEach(e => {
    if (e.isIntersecting) { e.target.classList.add("in"); io.unobserve(e.target); }
  }), { threshold: 0.12 });
  els.forEach(e => io.observe(e));
})();

// ---------------------------------------------------------------- small dense linear algebra
function zeros(n, m) { return Array.from({ length: n }, () => new Float64Array(m)); }
function luFactor(A) {
  const n = A.length, a = A.map(r => Float64Array.from(r)), p = new Int32Array(n);
  for (let k = 0; k < n; k++) {
    let mx = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(a[i][k]) > Math.abs(a[mx][k])) mx = i;
    [a[k], a[mx]] = [a[mx], a[k]]; p[k] = mx;
    for (let i = k + 1; i < n; i++) {
      const f = (a[i][k] /= a[k][k]);
      if (f) for (let j = k + 1; j < n; j++) a[i][j] -= f * a[k][j];
    }
  }
  return { a, p, n };
}
function luSolve(F, bIn) {
  const { a, p, n } = F, b = Float64Array.from(bIn);
  for (let k = 0; k < n; k++) { const t = b[k]; b[k] = b[p[k]]; b[p[k]] = t; }
  // the row swaps were applied to whole rows during the factorisation, so replay them on b
  for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) b[i] -= a[i][j] * b[j];
  for (let i = n - 1; i >= 0; i--) { for (let j = i + 1; j < n; j++) b[i] -= a[i][j] * b[j]; b[i] /= a[i][i]; }
  return b;
}
function matVec(A, x) { const n = A.length, y = new Float64Array(n); for (let i = 0; i < n; i++) { let s = 0; const r = A[i]; for (let j = 0; j < n; j++) s += r[j] * x[j]; y[i] = s; } return y; }

// ---------------------------------------------------------------- the stupa model
// Outline (height, half width) after the 70 m Phra Chedi Chai Mongkhon: two terraces, the
// ringed base, the bell, the throne and the spire. Vertical steps are repeated heights.
const PROFILE = [[0, 21], [4.3, 21], [4.3, 16], [12.8, 16], [12.8, 10.5], [15, 10.2], [17, 9.3], [19.5, 9.0], [21.5, 8.4],
  [23.5, 8.0], [26, 8.3], [30, 8.1], [34, 7.2], [37, 5.8], [39.5, 3.6], [42.5, 3.5], [42.5, 2.4], [47, 2.2], [52, 1.6],
  [56, 1.2], [60, 0.8], [64, 0.5], [69.7, 0.05]];
const H = 69.7;
function halfWidth(z) {
  for (let i = PROFILE.length - 1; i > 0; i--) {
    const [z0, r0] = PROFILE[i - 1], [z1, r1] = PROFILE[i];
    if (z >= z0 && z <= z1 && z1 > z0) return r0 + (r1 - r0) * (z - z0) / (z1 - z0);
  }
  return PROFILE[PROFILE.length - 1][1];
}

// A tapered cantilever of Hermite beam elements (v, θ per node), solid circular section of
// the outline's half width, masonry E 1.08 GPa (the modulus that matches the measured
// frequency of the real stupa in the calibration) and 1.8 t/m³. Units kN, m, t, s.
const NE = 24, LE = H / NE, E = 1.08e6, RHO = 1.8;
const ND = 2 * NE;                       // the base node is fixed
const K = zeros(ND, ND), M = zeros(ND, ND);
for (let e = 0; e < NE; e++) {
  const r = halfWidth((e + 0.5) * LE), A = Math.PI * r * r, I = Math.PI * r ** 4 / 4, L = LE;
  const k = [[12, 6 * L, -12, 6 * L], [6 * L, 4 * L * L, -6 * L, 2 * L * L], [-12, -6 * L, 12, -6 * L], [6 * L, 2 * L * L, -6 * L, 4 * L * L]];
  const m = [[156, 22 * L, 54, -13 * L], [22 * L, 4 * L * L, 13 * L, -3 * L * L], [54, 13 * L, 156, -22 * L], [-13 * L, -3 * L * L, -22 * L, 4 * L * L]];
  const ks = E * I / L ** 3, ms = RHO * A * L / 420;
  const dof = [2 * e - 2, 2 * e - 1, 2 * e, 2 * e + 1];
  for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
    if (dof[i] < 0 || dof[j] < 0) continue;
    K[dof[i]][dof[j]] += ks * k[i][j]; M[dof[i]][dof[j]] += ms * m[i][j];
  }
}
// the first frequency by inverse iteration
const KF = luFactor(K);
const f1 = (() => {
  let x = new Float64Array(ND).fill(1), lam = 0;
  for (let it = 0; it < 60; it++) {
    const y = luSolve(KF, matVec(M, x));
    const ky = matVec(K, y), my = matVec(M, y);
    let a = 0, b = 0; for (let i = 0; i < ND; i++) { a += y[i] * ky[i]; b += y[i] * my[i]; }
    lam = a / b;
    const n = Math.hypot(...y); x = y.map(v => v / n);
  }
  return Math.sqrt(lam) / (2 * Math.PI);
})();
// Rayleigh damping, 5 % at f1 and at the cantilever's third mode (≈ 6.27 f1 for a prism)
const w1 = 2 * Math.PI * f1, w2 = 6.27 * w1, XI = 0.05;
const ALPHA = 2 * XI * w1 * w2 / (w1 + w2), BETA = 2 * XI / (w1 + w2);
const C = K.map((row, i) => row.map((v, j) => ALPHA * M[i][j] + BETA * v));
const iota = new Float64Array(ND).map((_, i) => (i % 2 === 0 ? 1 : 0));
const Mi = matVec(M, iota);

// ground motions (m/s²), deterministic
function rng(seed) { return () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296); }
const DT = 0.01, DUR = 24;
const G = 9.81;
const MOTIONS = [
  { name: "Far-field motion", note: "long-period basin motion", pga: 0.05, make: (t, r) => {
    const env = t < 3 ? (t / 3) ** 2 : Math.exp(-(t - 3) / 7);
    let s = 0; for (const c of r) s += c.a * Math.sin(2 * Math.PI * c.f * t + c.p);
    return env * s;
  }, comps: seed => { const g = rng(seed); return Array.from({ length: 14 }, () => ({ f: 0.35 + 1.4 * g(), p: 6.283 * g(), a: 0.5 + g() })); } },
  { name: "Near-fault pulse", note: "a velocity pulse with high-frequency noise", pga: 0.22, make: (t, r) => {
    const tp = 1 / Math.max(0.6, f1 * 0.9), t0 = 4;
    const pulse = Math.sin(2 * Math.PI * (t - t0) / tp) * Math.exp(-(((t - t0) / (0.8 * tp)) ** 2));
    let n = 0; for (const c of r) n += c.a * Math.sin(2 * Math.PI * c.f * t + c.p);
    return pulse + 0.18 * n * Math.exp(-(((t - t0) / 2.5) ** 2));
  }, comps: seed => { const g = rng(seed); return Array.from({ length: 10 }, () => ({ f: 2 + 8 * g(), p: 6.283 * g(), a: g() })); } },
  { name: "Harmonic at f1", note: "resonance, then free vibration", pga: 0.01, make: t => {
    const ramp = t < 3 ? t / 3 : 1, on = t < 11 ? 1 : 0;
    return ramp * on * Math.sin(2 * Math.PI * f1 * t);
  }, comps: () => [] },
];

function response(mo) {
  const n = Math.round(DUR / DT), comps = mo.comps(7 + MOTIONS.indexOf(mo) * 31);
  const ag = new Float64Array(n + 1);
  for (let i = 0; i <= n; i++) ag[i] = mo.make(i * DT, comps);
  const peak = Math.max(...ag.map(Math.abs)); for (let i = 0; i <= n; i++) ag[i] *= mo.pga * G / peak;
  // Newmark, average acceleration
  const a0 = 4 / (DT * DT), a1 = 2 / DT;
  const Keff = K.map((row, i) => row.map((v, j) => v + a0 * M[i][j] + a1 * C[i][j]));
  const F = luFactor(Keff);
  let u = new Float64Array(ND), v = new Float64Array(ND), a = new Float64Array(ND);
  const frames = [u];
  for (let s = 1; s <= n; s++) {
    const m1 = new Float64Array(ND), c1 = new Float64Array(ND);
    for (let i = 0; i < ND; i++) { m1[i] = a0 * u[i] + 4 / DT * v[i] + a[i]; c1[i] = a1 * u[i] + v[i]; }
    const Mm = matVec(M, m1), Cc = matVec(C, c1);
    const p = new Float64Array(ND); for (let i = 0; i < ND; i++) p[i] = -Mi[i] * ag[s] + Mm[i] + Cc[i];
    const un = luSolve(F, p), an = new Float64Array(ND), vn = new Float64Array(ND);
    for (let i = 0; i < ND; i++) { an[i] = a0 * (un[i] - u[i]) - 4 / DT * v[i] - a[i]; vn[i] = v[i] + DT / 2 * (a[i] + an[i]); }
    u = un; v = vn; a = an; frames.push(u);
  }
  return { ag, frames };
}

// displacement v, rotation θ and curvature κ at height z from the nodal values
function beamAt(u, z) {
  const e = Math.min(NE - 1, Math.max(0, Math.floor(z / LE))), x = z / LE - e, L = LE;
  const g = i => (i < 0 ? 0 : u[i]);
  const v1 = g(2 * e - 2), t1 = g(2 * e - 1), v2 = g(2 * e), t2 = g(2 * e + 1);
  const N1 = 1 - 3 * x * x + 2 * x ** 3, N2 = L * (x - 2 * x * x + x ** 3), N3 = 3 * x * x - 2 * x ** 3, N4 = L * (-x * x + x ** 3);
  const d1 = (-6 * x + 6 * x * x) / L, d2 = 1 - 4 * x + 3 * x * x, d3 = (6 * x - 6 * x * x) / L, d4 = -2 * x + 3 * x * x;
  const s1 = (-6 + 12 * x) / (L * L), s2 = (-4 + 6 * x) / L, s3 = (6 - 12 * x) / (L * L), s4 = (-2 + 6 * x) / L;
  return [N1 * v1 + N2 * t1 + N3 * v2 + N4 * t2, d1 * v1 + d2 * t1 + d3 * v2 + d4 * t2, s1 * v1 + s2 * t1 + s3 * v2 + s4 * t2];
}

// the meshed outline: bands along the profile, 12 cells across, two triangles per cell
const MESH = (() => {
  const verts = [], tris = [];
  const NX = 12;
  for (let i = 1; i < PROFILE.length; i++) {
    const [z0, r0] = PROFILE[i - 1], [z1, r1] = PROFILE[i];
    if (z1 <= z0) continue;
    const nb = Math.max(1, Math.round((z1 - z0) / 1.6));
    const base = verts.length;
    for (let b = 0; b <= nb; b++) {
      const z = z0 + (z1 - z0) * b / nb, r = r0 + (r1 - r0) * b / nb;
      for (let c = 0; c <= NX; c++) verts.push([(-1 + 2 * c / NX) * r, z]);
    }
    for (let b = 0; b < nb; b++) for (let c = 0; c < NX; c++) {
      const p = base + b * (NX + 1) + c, q = p + NX + 1;
      if ((b + c) % 2) { tris.push([p, p + 1, q + 1], [p, q + 1, q]); } else { tris.push([p, p + 1, q], [p + 1, q + 1, q]); }
    }
  }
  return { verts, tris };
})();

// the contour colours of the program's legend: 12 bands, blue to red
const BANDS = ["#0000ff", "#0046ff", "#008cff", "#00d2ff", "#00ffb4", "#00ff46", "#46ff00", "#b4ff00", "#ffd200", "#ff8c00", "#ff4600", "#ff0000"];
function bandColor(t) { return BANDS[Math.max(0, Math.min(11, Math.floor((t + 1) / 2 * 12)))]; }

// ---------------------------------------------------------------- hero animation
(() => {
  const cv = document.getElementById("sim"), tr = document.getElementById("trace");
  if (!cv) return;
  const ctx = cv.getContext("2d"), tctx = tr.getContext("2d");
  const stage = document.getElementById("simStage"), timeEl = document.getElementById("simTime");
  const driftEl = document.getElementById("simDrift"), freqEl = document.getElementById("simFreq");
  const legend = document.getElementById("legendBands");
  freqEl.textContent = `f1 ${f1.toFixed(2)} Hz`;
  let mi = 0, run = null, frame = 0, scale = 1, smax = 1;

  function fit(c) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2), w = c.clientWidth, h = c.clientHeight;
    if (c.width !== Math.round(w * dpr) || c.height !== Math.round(h * dpr)) { c.width = Math.round(w * dpr); c.height = Math.round(h * dpr); }
    return [c.width, c.height, dpr];
  }
  function start(i) {
    mi = i; const mo = MOTIONS[i]; run = response(mo); frame = 0;
    const top = 2 * NE - 2;
    const vmax = Math.max(...run.frames.map(u => Math.abs(u[top]))) || 1e-9;
    scale = 0.07 * H / vmax;
    smax = 1e-9;
    for (let s = 0; s < run.frames.length; s += 5) for (let z = 0.2; z < H; z += 1.5) {
      const k = beamAt(run.frames[s], z)[2]; smax = Math.max(smax, Math.abs(E * k * halfWidth(z)));
    }
    lim = smax * 0.05;
    stage.textContent = `${mo.name} · PGA ${mo.pga} g`;
  }
  // contour limits per frame, as the program's auto-scaled legend, eased so they do not flicker
  let lim = 1, legendAt = -1;
  function setLegend() {
    legend.innerHTML = BANDS.slice().reverse().map((c, k) => {
      const val = lim * (1 - 2 * k / 11);
      return `<div><i style="background:${c}"></i>${val >= 0 ? "+" : ""}${val.toExponential(2)}</div>`;
    }).join("");
  }
  document.getElementById("simNext").addEventListener("click", () => start((mi + 1) % MOTIONS.length));

  function draw() {
    const [W, Hh, dpr] = fit(cv);
    ctx.clearRect(0, 0, W, Hh);
    const u = run.frames[Math.min(frame, run.frames.length - 1)];
    const sc = Math.min((Hh * 0.74) / H, (W * 0.62) / (2 * 21 + 10));
    const ox = W * 0.56, oy = Hh * 0.83;
    // ground
    ctx.fillStyle = "#c9cfdf"; ctx.fillRect(0, oy, W, Hh - oy);
    ctx.strokeStyle = "#9aa3bb"; ctx.lineWidth = 1 * dpr;
    for (let x = -30; x < W / dpr + 30; x += 14) { ctx.beginPath(); ctx.moveTo(x * dpr, oy + 12 * dpr); ctx.lineTo((x + 10) * dpr, oy + 2 * dpr); ctx.stroke(); }
    // deformed, coloured mesh
    const P = MESH.verts.map(([x, z]) => {
      const [v, th] = beamAt(u, z);
      return [ox + (x + scale * v) * sc, oy - (z - scale * x * th) * sc];
    });
    ctx.lineWidth = 0.5 * dpr; ctx.strokeStyle = "rgba(20,24,51,0.35)";
    const sig = MESH.tris.map(([a, b, c]) => {
      const xc = (MESH.verts[a][0] + MESH.verts[b][0] + MESH.verts[c][0]) / 3;
      const zc = (MESH.verts[a][1] + MESH.verts[b][1] + MESH.verts[c][1]) / 3;
      return -E * beamAt(u, zc)[2] * xc;              // S33 = −E κ x (tension positive)
    });
    const now = Math.max(smax * 0.02, ...sig.map(Math.abs));
    lim = now > lim ? now : lim + 0.06 * (now - lim);
    if (frame - legendAt >= 10 || frame < legendAt) { legendAt = frame; setLegend(); }
    for (let i = 0; i < MESH.tris.length; i++) {
      const [a, b, c] = MESH.tris[i];
      ctx.fillStyle = bandColor(sig[i] / lim);
      ctx.beginPath(); ctx.moveTo(...P[a]); ctx.lineTo(...P[b]); ctx.lineTo(...P[c]); ctx.closePath(); ctx.fill(); ctx.stroke();
    }
    // labels like the program's view
    ctx.fillStyle = "#141833"; ctx.font = `${11 * dpr}px "JetBrains Mono", monospace`;
    const t = frame * DT;
    ctx.fillStyle = "rgba(255,255,255,0.78)"; ctx.fillRect(8 * dpr, Hh - 54 * dpr, 430 * dpr, 46 * dpr);
    ctx.fillStyle = "#141833";
    ctx.fillText(`Step: Earthquake   t = ${t.toFixed(2)} s`, 14 * dpr, Hh - 40 * dpr);
    ctx.fillText(`Primary Var: S, S33`, 14 * dpr, Hh - 27 * dpr);
    ctx.fillText(`Deformed Var: U   Deformation Scale Factor: +${scale.toExponential(2)}`, 14 * dpr, Hh - 14 * dpr);
    timeEl.textContent = `t = ${t.toFixed(1)} s`;
    driftEl.textContent = `top drift ${(1000 * u[2 * NE - 2]).toFixed(1)} mm`;

    // the ground record with a cursor
    const [tw, th2] = fit(tr);
    tctx.clearRect(0, 0, tw, th2);
    const ag = run.ag, n = ag.length, amax = MOTIONS[mi].pga * G;
    tctx.strokeStyle = "#c3c8dc"; tctx.beginPath(); tctx.moveTo(0, th2 / 2); tctx.lineTo(tw, th2 / 2); tctx.stroke();
    tctx.strokeStyle = "#3a3fd6"; tctx.lineWidth = 1.2 * dpr; tctx.beginPath();
    for (let i = 0; i < n; i += 2) { const x = i / n * tw, y = th2 / 2 - ag[i] / amax * th2 * 0.42; i ? tctx.lineTo(x, y) : tctx.moveTo(x, y); }
    tctx.stroke();
    const cx = frame / n * tw;
    tctx.fillStyle = "rgba(58,63,214,0.10)"; tctx.fillRect(0, 0, cx, th2);
    tctx.strokeStyle = "#f2a51a"; tctx.lineWidth = 2 * dpr; tctx.beginPath(); tctx.moveTo(cx, 0); tctx.lineTo(cx, th2); tctx.stroke();
    tctx.fillStyle = "#4a5070"; tctx.font = `${10 * dpr}px "JetBrains Mono", monospace`;
    tctx.fillText(`ground acceleration · ${MOTIONS[mi].note}`, 6 * dpr, 12 * dpr);
  }

  start(0);
  let last = performance.now(), acc = 0;
  function loop(now) {
    acc += Math.min(0.1, (now - last) / 1000); last = now;
    while (acc >= DT) { acc -= DT; frame++; }
    if (frame >= run.frames.length) start((mi + 1) % MOTIONS.length);
    draw();
    requestAnimationFrame(loop);
  }
  if (reduceMotion) { frame = Math.round(6 / DT); draw(); } else requestAnimationFrame(loop);
  window.addEventListener("resize", draw);
})();

// ---------------------------------------------------------------- replayed agent session
const TOOLS = ["current", "new_model", "summary", "get_model", "set_model", "add_material", "add_section", "add_part", "select", "add_set",
  "add_bc", "add_load", "add_step", "add_amplitude", "add_tie", "add_spring", "add_mass", "add_monitor", "delete", "mesh", "check", "run",
  "cancel", "get_results", "get_history", "show", "output", "screenshot", "set_compute", "get_compute", "workdir", "save", "open", "export"];

const SCRIPT = [
  { user: "Model a 4 m pinned JIS H-200×200×8×12 steel column under 1,000 kN and find its buckling load. Show me the first mode." },
  { say: "I'll check what the window has open first." },
  { tool: "current", args: "", out: "No file open: an empty model in the window.", view: "empty" },
  { tool: "new_model", args: 'name: "H-200 column, pinned", replace: true', out: "New model 'H-200 column, pinned'." },
  { tool: "add_material", args: 'name: "SS400", E: 2.05e8, nu: 0.3, density: 7.85', out: "Material SS400 (elastic)." },
  { tool: "add_section", args: 'name: "H200", material: "SS400", type: "beam", profile: {shape: "i", h: 0.2, b: 0.2, tw: 0.008, tf: 0.012}', out: "Beam section H200." },
  { tool: "add_part", args: 'name: "column", section: "H200", shape: {type: "line", points: [[0,0,0],[0,0,4]]}, mesh_size: 0.25', out: "Part column: 16 beam elements.", view: "column" },
  { tool: "add_bc", args: 'name: "base", at: {near: [[0,0,0]]}, dofs: "1,2,3,6"', out: "Support base on 1 node.", view: "base" },
  { tool: "add_bc", args: 'name: "top", at: {near: [[0,0,4]]}, dofs: "1,2"', out: "Support top on 1 node.", view: "top" },
  { tool: "add_step", args: 'name: "Buckle", type: "buckle", modes: 3', out: "Step Buckle (linear buckling, 3 modes)." },
  { tool: "add_load", args: 'name: "P", type: "force", at: {near: [[0,0,4]]}, f: [0,0,-1000], step: "Buckle"', out: "Load P: the reference load of Buckle.", view: "load" },
  { tool: "run", args: "", out: "Buckle (buckle): completed, 3 frames. Buckling load factor 1 = 2.0113", view: "run" },
  { tool: "show", args: 'module: "output", field: "U", step: "Buckle", frame: 0', out: "The window shows U, Magnitude, buckling mode 1.", view: "mode" },
  { say: "The column buckles about its weak axis at **λ = 2.011**, a critical load of **2,011 kN**. Euler's π²EI/L² with I = 1.60×10⁻⁵ m⁴ gives 2,024 kN, so the beam model is 0.6 % low. Mode 1 is on screen; the model is unsaved, say where to save it." },
];

(() => {
  const body = document.getElementById("termBody"), grid = document.getElementById("toolGrid");
  const mini = document.getElementById("mini"), tag = document.getElementById("viewTag");
  if (!body) return;
  grid.innerHTML = TOOLS.map(t => `<span data-t="${t}">${t}</span>`).join("");
  const mctx = mini.getContext("2d");
  let state = "empty", animT = 0, timers = [], playing = false;

  function drawMini() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2), w = mini.clientWidth, h = mini.clientHeight;
    mini.width = w * dpr; mini.height = h * dpr;
    const c = mctx; c.clearRect(0, 0, mini.width, mini.height); c.save(); c.scale(dpr, dpr);
    const order = ["empty", "column", "base", "top", "load", "run", "mode"], lvl = order.indexOf(state);
    const x0 = w / 2, yb = h - 34, yt = 40, L = yb - yt;
    c.fillStyle = "#4a5070"; c.font = '11px "JetBrains Mono", monospace';
    if (lvl < 1) { c.fillText("empty model", x0 - 38, h / 2); c.restore(); return; }
    const amp = lvl >= 6 ? 26 * Math.sin(animT * 2.2) : 0;
    const N = 16;
    const pts = Array.from({ length: N + 1 }, (_, i) => { const s = i / N; return [x0 + amp * Math.sin(Math.PI * s), yb - s * L, Math.abs(Math.sin(Math.PI * s))]; });
    if (lvl >= 6) {
      for (let i = 0; i < N; i++) {
        c.strokeStyle = bandColor(2 * (pts[i][2] + pts[i + 1][2]) / 2 * Math.abs(Math.sin(animT * 2.2)) - 1); c.lineWidth = 9;
        c.beginPath(); c.moveTo(pts[i][0], pts[i][1]); c.lineTo(pts[i + 1][0], pts[i + 1][1]); c.stroke();
      }
    } else {
      c.strokeStyle = "#8f96b3"; c.lineWidth = 9; c.beginPath(); c.moveTo(x0, yb); c.lineTo(x0, yt); c.stroke();
      c.fillStyle = "#141833"; for (const p of pts) { c.beginPath(); c.arc(p[0], p[1], 2.2, 0, 7); c.fill(); }
    }
    const pin = (y, up) => {
      c.fillStyle = "#f2a51a"; c.beginPath();
      if (up) { c.moveTo(x0, y); c.lineTo(x0 - 9, y + 14); c.lineTo(x0 + 9, y + 14); } else { c.moveTo(x0, y); c.lineTo(x0 - 9, y - 14); c.lineTo(x0 + 9, y - 14); }
      c.closePath(); c.fill();
    };
    if (lvl >= 2) pin(yb + 2, true);
    if (lvl >= 3) { c.fillStyle = "#f2a51a"; c.fillRect(x0 - 22, yt - 2, 10, 4); c.fillRect(x0 + 12, yt - 2, 10, 4); }
    if (lvl >= 4) {
      c.strokeStyle = "#e5484d"; c.fillStyle = "#e5484d"; c.lineWidth = 2.5;
      c.beginPath(); c.moveTo(x0, yt - 34); c.lineTo(x0, yt - 6); c.stroke();
      c.beginPath(); c.moveTo(x0, yt - 2); c.lineTo(x0 - 6, yt - 12); c.lineTo(x0 + 6, yt - 12); c.closePath(); c.fill();
      c.fillStyle = "#141833"; c.fillText("1000 kN", x0 + 10, yt - 22);
    }
    if (lvl === 5) { c.fillStyle = "#3a3fd6"; c.fillText("solving…", x0 + 18, h / 2); }
    if (lvl >= 6) { c.fillStyle = "#141833"; c.fillText("Buckling mode 1: λ = 2.0113", 12, 18); }
    c.restore();
  }
  const tags = { empty: "empty model", column: "H200 column · 16 B31", base: "support: base", top: "support: top", load: "load P = 1000 kN", run: "running Buckle", mode: "Output · U, Magnitude · mode 1" };
  function setView(v) { state = v; tag.textContent = tags[v]; drawMini(); }
  (function spin() { animT += 0.016; if (state === "mode") drawMini(); requestAnimationFrame(spin); })();

  const esc = s => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/\*\*(.+?)\*\*/g, "<b>$1</b>");
  function line(cls, html) { const d = document.createElement("div"); d.className = cls; d.innerHTML = html; body.appendChild(d); body.scrollTop = body.scrollHeight; return d; }
  function later(ms, f) { timers.push(setTimeout(f, ms)); }
  function play() {
    timers.forEach(clearTimeout); timers = []; body.innerHTML = ""; playing = true;
    grid.querySelectorAll("span").forEach(s => s.classList.remove("hit", "now"));
    setView("empty");
    line("t-head", "✻ Claude Code · MCP server continuax connected to the window");
    let t = 500;
    for (const st of SCRIPT) {
      if (st.user) {
        const d = line("t-user", ""); const txt = st.user;
        for (let i = 0; i <= txt.length; i += 2) later(t + i * 14, () => { d.textContent = txt.slice(0, i); });
        later(t + txt.length * 14 + 20, () => { d.textContent = txt; });
        t += txt.length * 14 + 700;
      } else if (st.say) {
        later(t, () => line("t-say", "⏺ " + esc(st.say))); t += 900;
      } else {
        later(t, () => {
          line("t-dot", `<b>⏺</b> <span class="t-tool">continuax - ${st.tool} (MCP)</span>(${esc(st.args)})`);
          grid.querySelectorAll("span").forEach(s => s.classList.remove("now"));
          const chip = grid.querySelector(`[data-t="${st.tool}"]`); if (chip) chip.classList.add("hit", "now");
        });
        t += st.tool === "run" ? 1500 : 650;
        later(t, () => { line("t-out", "  ⎿  " + esc(st.out)); if (st.view) setView(st.view); });
        t += 450;
      }
    }
    later(t, () => { line("", '<span class="cursor"></span>'); playing = false; });
  }
  document.getElementById("replay").addEventListener("click", play);
  const io = "IntersectionObserver" in window ? new IntersectionObserver(es => {
    if (es.some(e => e.isIntersecting) && !playing && !body.childElementCount) play();
  }, { threshold: 0.3 }) : null;
  if (io) io.observe(document.getElementById("term")); else play();
  window.addEventListener("resize", drawMini);
  drawMini();
})();

// trial and quotation request: builds an email to tellusfem@gmail.com (the site has no server)
(function () {
  const form = document.getElementById("reqForm");
  if (!form) return;
  const TO = "tellusfem@gmail.com";
  const send = document.getElementById("reqSend"), err = document.getElementById("reqErr"), note = document.getElementById("reqNote");
  const val = (n) => (form.elements[n] ? String(form.elements[n].value || "").trim() : "");
  const isQuote = () => form.elements.type.value === "quote";
  function sync() {
    form.classList.toggle("is-quote", isQuote());
    send.textContent = isQuote() ? "Request a quotation" : "Request the 30-day trial";
  }
  form.addEventListener("change", sync);
  sync();
  function build() {
    const kind = isQuote() ? "Quotation request" : "30-day trial request";
    const rows = [
      ["Request", isQuote() ? "Quotation" : "30-day trial"],
      ["Name", val("name")], ["Email", val("email")], ["Organisation", val("org")], ["Country", val("country")],
      ["Main use", val("use")], ["Computer", val("platform")],
    ];
    if (isQuote()) rows.push(["Number of users", val("seats")], ["Licence", val("licence")]);
    const body = rows.filter((r) => r[1]).map((r) => r[0] + ": " + r[1]).join("\n") + (val("msg") ? "\n\nMessage:\n" + val("msg") : "");
    return { subject: "ContinuaX " + kind + " - " + (val("org") || val("name")), body };
  }
  function check() {
    let bad = null;
    ["name", "email", "org"].forEach((n) => {
      const el = form.elements[n], ok = el.value.trim() !== "" && (n !== "email" || /^[^@\s]+@[^@\s]+\.[^@\s]+$/.test(el.value.trim()));
      el.setAttribute("aria-invalid", ok ? "false" : "true");
      if (!ok && !bad) bad = el;
    });
    if (bad) { err.textContent = "Please fill in your name, a valid work email and your organisation."; bad.focus(); return false; }
    err.textContent = ""; return true;
  }
  form.addEventListener("submit", (e) => {
    e.preventDefault();
    if (!check()) return;
    const m = build();
    window.location.href = "mailto:" + TO + "?subject=" + encodeURIComponent(m.subject) + "&body=" + encodeURIComponent(m.body);
    note.textContent = "Your mail app should now show the request to " + TO + ". Press send there to finish.";
  });
  document.getElementById("reqCopy").addEventListener("click", async () => {
    if (!check()) return;
    const m = build(), text = "To: " + TO + "\nSubject: " + m.subject + "\n\n" + m.body;
    try { await navigator.clipboard.writeText(text); note.textContent = "Copied. Paste it into an email to " + TO + "."; }
    catch (_) { note.textContent = "Copy is blocked in this browser. Please write to " + TO + "."; }
  });
})();
