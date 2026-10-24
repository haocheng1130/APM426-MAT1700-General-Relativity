// Integrate the example orbits of kerr_geodesics.html with the page's own physics code and print them as JSON.
// Run with Node.js: node dump_trajectories.js > trajectories.json
const fs = require('fs');
const path = require('path');

const html = fs.readFileSync(path.join(__dirname, '..', 'kerr_geodesics.html'), 'utf8');
const block = html.match(/\/\* ===== physics: begin ===== \*\/([\s\S]*?)\/\* ===== physics: end ===== \*\//);
const Kerr = new Function(block[1] + '\nreturn Kerr;')();
const deg = Math.PI / 180;
const sp = Kerr.sphericalPhoton(0.9, 2.6);

// The same examples as the page's PRESETS list.
const examples = [
  { name: 'inclined', a: 0.9, rp: 6, ra: 12, inc: 45, tMax: 4000, expect: 'done' },
  { name: 'polar', a: 0.9, rp: 7, ra: 11, inc: 90, tMax: 4000, expect: 'done' },
  { name: 'retrograde', a: 0.9, rp: 9, ra: 14, inc: 155, tMax: 4000, expect: 'done' },
  { name: 'zoom-whirl', a: 0.7, rp: 2.68, ra: 16, inc: 20, tMax: 2500, expect: 'done' },
  { name: 'schwarzschild', a: 0, rp: 8, ra: 20, inc: 0, tMax: 6000, expect: 'done' },
  { name: 'plunge', a: 0.9, k: -1, E: 0.97, L: 2.2, Q: 3, r0: 14, th0: 80, sr: -1, sth: 1, tMax: 1000, expect: 'captured' },
  { name: 'photon sphere', a: 0.9, k: 0, E: 1, L: sp.L, Q: sp.Q, r0: 2.6, th0: 90, sr: 1, sth: 1, tMax: 400, expect: 'escaped' },
  { name: 'light bending', a: 0.9, k: 0, E: 1, L: 5, Q: 12, r0: 40, th0: 90, sr: -1, sth: 1, tMax: 400, expect: 'escaped' },
];

const out = examples.map(ex => {
  let c, y0, rEscape;
  if (ex.rp !== undefined) {
    c = Kerr.solveBound(ex.a, ex.rp, ex.ra, Math.cos(ex.inc * deg));
    y0 = [ex.rp, 0, Math.PI / 2, Math.sqrt(Math.max(Kerr.Theta(c, Math.PI / 2), 0)), 0, 0, 0];
    rEscape = 4 * ex.ra + 20;
  } else {
    c = { a: ex.a, E: ex.E, L: ex.L, Q: ex.Q, k: ex.k };
    const th0 = ex.th0 * deg;
    y0 = [ex.r0, ex.sr * Math.sqrt(Math.max(Kerr.R(c, ex.r0), 0)), th0, ex.sth * Math.sqrt(Math.max(Kerr.Theta(c, th0), 0)), 0, 0, 0];
    rEscape = Math.max(80, 2.5 * ex.r0);
  }
  const tr = Kerr.integrate(c, y0, { tMax: ex.tMax, rEscape });
  const n = tr.t.length;
  const samples = Array.from({ length: 60 }, (_, i) => Math.floor(i * (n - 1) / 59))
    .map(i => [tr.lam[i], tr.x[i], tr.y[i], tr.z[i], tr.t[i]]);
  return { name: ex.name, expect: ex.expect, status: tr.status, drift: tr.drift, a: c.a, E: c.E, L: c.L, y0, samples };
});
process.stdout.write(JSON.stringify(out));
