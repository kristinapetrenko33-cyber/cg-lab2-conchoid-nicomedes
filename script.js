'use strict';
/* ======================================================================
 *  Лабораторна робота №2 · Побудова чудових кривих · Варіант 10
 *  Конхоїда Нікомеда   (x − a)²(x² + y²) = b²x²       Петренко К. Ю.
 *
 *  Запис кривої — параметричний (формула (7) методички), параметр φ:
 *        x = a + b·cos φ,    y = a·tg φ + b·sin φ
 *  (еквівалент полярного r = a / cos φ + b). Похідні виведено вручну:
 *        x′ = −b sin φ,                 y′ = a/cos²φ + b cos φ
 *        x″ = −b cos φ,                 y″ = 2a sin φ/cos³φ − b sin φ
 *  Чисельно: адаптивний метод Сімпсона (довжини, площі), многочлен
 *  b·c³ + 3a·c² − 2a = 0, c = cos φ, — для точок перегину.
 *  Рендер: WebGL (власна тесселяція ліній) або резервний Canvas 2D.
 * ==================================================================== */

const PI = Math.PI, TAU = 2 * PI;
const $ = id => document.getElementById(id);
const rad = d => d * PI / 180, deg = r => r * 180 / PI;
const num = id => { const v = parseFloat($(id).value); return Number.isFinite(v) ? v : 0; };

/* ---------------------------------------------------- 1. математика */
const K = {
    pt: (p, a, b) => [a + b * Math.cos(p), a * Math.tan(p) + b * Math.sin(p)],
    d1: (p, a, b) => { const c = Math.cos(p); return [-b * Math.sin(p), a / (c * c) + b * c]; },
    d2: (p, a, b) => { const c = Math.cos(p); return [-b * c, 2 * a * Math.sin(p) / (c * c * c) - b * Math.sin(p)]; },
    speed(p, a, b) { return Math.hypot(...K.d1(p, a, b)); },
    cross(p, a, b) { const [x1, y1] = K.d1(p, a, b), [x2, y2] = K.d2(p, a, b); return x1 * y2 - y1 * x2; },
    curv(p, a, b) { return K.cross(p, a, b) / K.speed(p, a, b) ** 3; },
    r: (p, a, b) => a / Math.cos(p) + b
};
const BRANCHES = [[-PI / 2, PI / 2], [PI / 2, 3 * PI / 2]];

/** адаптивний Сімпсон із оцінкою похибки за Річардсоном */
function adaptiveSimpson(f, lo, hi, eps = 1e-9) {
    const rec = (a, b, e, whole, fa, fm, fb, depth) => {
        const m = (a + b) / 2, l = (a + m) / 2, r = (m + b) / 2, fl = f(l), fr = f(r);
        const left = (m - a) / 6 * (fa + 4 * fl + fm), right = (b - m) / 6 * (fm + 4 * fr + fb);
        if (depth <= 0 || Math.abs(left + right - whole) <= 15 * e) return left + right + (left + right - whole) / 15;
        return rec(a, m, e / 2, left, fa, fl, fm, depth - 1) + rec(m, b, e / 2, right, fm, fr, fb, depth - 1);
    };
    const m = (lo + hi) / 2, fa = f(lo), fm = f(m), fb = f(hi);
    return rec(lo, hi, eps, (hi - lo) / 6 * (fa + 4 * fm + fb), fa, fm, fb, 22);
}
/** чи міститься полюс cos φ = 0 у відкритому проміжку (lo, hi) */
const crossesPole = (lo, hi) => {
    for (let k = Math.ceil((lo - PI / 2) / PI); PI / 2 + k * PI < hi; k++) if (PI / 2 + k * PI > lo) return true;
    return false;
};
const arcLength = (lo, hi, a, b) => crossesPole(lo, hi) ? Infinity : adaptiveSimpson(p => K.speed(p, a, b), lo, hi);
const loopRange = (a, b) => b > a ? (f0 => [PI - f0, PI + f0])(Math.acos(a / b)) : null;
const loopAreaNumeric = (a, b) => { const r = loopRange(a, b); return r ? 0.5 * adaptiveSimpson(p => K.r(p, a, b) ** 2, r[0], r[1]) : 0; };
const loopAreaExact = (a, b) => {
    if (b <= a) return 0;
    const f0 = Math.acos(a / b);
    return b * b * f0 - 2 * a * b * Math.log(Math.tan(PI / 4 + f0 / 2)) + a * a * Math.tan(f0);
};
/** корені g(c) = b c³ + 3a c² − 2a на (−1; 1) → відповідні φ ∈ [−π/2; 3π/2) */
function inflections(a, b) {
    const g = c => b * c ** 3 + 3 * a * c * c - 2 * a, out = [];
    let c0 = -1 + 1e-9, g0 = g(c0);
    for (let i = 1; i <= 4000; i++) {
        const c1 = -1 + 2 * i / 4000, g1 = g(c1);
        if (g0 * g1 < 0) {
            let lo = c0, hi = c1;
            for (let k = 0; k < 60; k++) { const m = (lo + hi) / 2; (g(lo) * g(m) <= 0) ? hi = m : lo = m; }
            const c = (lo + hi) / 2, ac = Math.acos(c);
            for (const p of c > 0 ? [ac, -ac] : [ac, TAU - ac]) out.push(p);
        }
        c0 = c1; g0 = g1;
    }
    return out.sort((u, v) => u - v);
}

/* ---------------------------------- 2. матриці (column-major) */
const Mx = {
    mul(A, B) { const o = new Array(9); for (let c = 0; c < 3; c++) for (let r = 0; r < 3; r++) o[c * 3 + r] = A[r] * B[c * 3] + A[3 + r] * B[c * 3 + 1] + A[6 + r] * B[c * 3 + 2]; return o; },
    translate: (dx, dy) => [1, 0, 0, 0, 1, 0, dx, dy, 1],
    rotate(d, px, py) { const c = Math.cos(rad(d)), s = Math.sin(rad(d)); return [c, s, 0, -s, c, 0, px - c * px + s * py, py - s * px - c * py, 1]; },
    apply: (m, [x, y]) => [m[0] * x + m[3] * y + m[6], m[1] * x + m[4] * y + m[7]]
};

/* ---------------------------------------------------------- 3. стан */
const S = {
    a: 2, b: 3, phi: 40,
    shift: [0, 0], rot: { x: 0, y: 0, a: 0 },
    view: { ox: 0, oy: 0, sc: 42 },
    show: { tan: true, nor: true, asy: true, circ: false, loop: true },
    calc: [-60, 60],
    anim: { on: false, t0: 0, speed: 4 },
    samples: []
};
const eucl = () => Mx.mul(Mx.rotate(S.rot.a, S.rot.x, S.rot.y), Mx.translate(S.shift[0], S.shift[1]));
const toPx = w => [S.view.ox + w[0] * S.view.sc, S.view.oy + w[1] * S.view.sc];
const toWorld = (x, y) => [(x - S.view.ox) / S.view.sc, (y - S.view.oy) / S.view.sc];

/* ---------------------------------------------------- 4. рендерер */
const STYLE = {
    grid:  { c: [0.15, 0.33, 0.55, 1], w: 1 },
    major: { c: [0.22, 0.44, 0.70, 1], w: 1 },
    axisX: { c: [1.0, 0.42, 0.51, 1], w: 1.8 },
    axisY: { c: [0.49, 0.91, 0.53, 1], w: 1.8 },
    curve: { c: [1, 1, 1, 1], w: 2.8 },
    asy:   { c: [1, 0.55, 0.55, 1], w: 1.5, dash: [9, 6] },
    tan:   { c: [0.5, 0.95, 0.6, 1], w: 1.8 },
    nor:   { c: [0.45, 0.82, 1, 1], w: 1.8 },
    circ:  { c: [1, 0.82, 0.4, 1], w: 1.3, dash: [5, 5] },
    pivot: { c: [1, 0.42, 0.5, 1], w: 1.6 },
    loop:  { c: [0.45, 0.82, 1, 0.20], w: 0 },
    dot:   { c: [1, 0.82, 0.4, 1], w: 0 }
};
const ORDER = ['loop', 'grid', 'major', 'axisX', 'axisY', 'asy', 'circ', 'curve', 'tan', 'nor', 'pivot', 'dot'];
const FILLS = new Set(['loop', 'dot']);
const glc = $('gl'), ov = $('ov'), g2 = ov.getContext('2d');
let gl = null, prog, buf, loc, W = 900, H = 600, dpr = 1;

function initGL() {
    try { gl = glc.getContext('webgl', { antialias: true, alpha: false }); } catch (e) { gl = null; }
    if (!gl) return;
    const sh = (t, s) => { const o = gl.createShader(t); gl.shaderSource(o, s); gl.compileShader(o); return o; };
    prog = gl.createProgram();
    gl.attachShader(prog, sh(gl.VERTEX_SHADER, 'attribute vec2 p; uniform vec2 res; void main(){ gl_Position = vec4(p / res * 2.0 - 1.0, 0.0, 1.0); }'));
    gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, 'precision mediump float; uniform vec4 col; void main(){ gl_FragColor = col; }'));
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) { gl = null; return; }
    gl.useProgram(prog); buf = gl.createBuffer();
    loc = { p: gl.getAttribLocation(prog, 'p'), res: gl.getUniformLocation(prog, 'res'), col: gl.getUniformLocation(prog, 'col') };
    gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
}

/** відсікання відрізка прямокутником екрана (Ліанга–Барські) */
function clipSeg(p, q, m = 40) {
    let t0 = 0, t1 = 1; const dx = q[0] - p[0], dy = q[1] - p[1];
    for (const [pp, qq] of [[-dx, p[0] + m], [dx, W + m - p[0]], [-dy, p[1] + m], [dy, H + m - p[1]]]) {
        if (pp === 0) { if (qq < 0) return null; }
        else { const r = qq / pp; if (pp < 0) { if (r > t1) return null; if (r > t0) t0 = r; } else { if (r < t0) return null; if (r < t1) t1 = r; } }
    }
    return [[p[0] + t0 * dx, p[1] + t0 * dy], [p[0] + t1 * dx, p[1] + t1 * dy]];
}
/** ламана → набір неперервних частин, що потрапили в кадр */
function runs(pts) {
    const out = []; let cur = null;
    for (let i = 1; i < pts.length; i++) {
        const c = clipSeg(pts[i - 1], pts[i]);
        if (!c) { cur = null; continue; }
        if (cur && Math.hypot(cur[cur.length - 1][0] - c[0][0], cur[cur.length - 1][1] - c[0][1]) < 1e-6) cur.push(c[1]);
        else { cur = [c[0], c[1]]; out.push(cur); }
    }
    return out;
}
function dashed(pts, pat) {
    const out = []; let cur = [pts[0]], i = 0, left = pat[0], on = true;
    for (let k = 1; k < pts.length; k++) {
        let [x0, y0] = pts[k - 1]; const [x1, y1] = pts[k];
        let len = Math.hypot(x1 - x0, y1 - y0); const ux = (x1 - x0) / (len || 1), uy = (y1 - y0) / (len || 1);
        while (len > left) {
            x0 += ux * left; y0 += uy * left; len -= left;
            if (on) { cur.push([x0, y0]); out.push(cur); cur = []; } else cur = [[x0, y0]];
            on = !on; i = (i + 1) % pat.length; left = pat[i];
        }
        left -= len; if (on) cur.push([x1, y1]);
    }
    if (on && cur.length > 1) out.push(cur);
    return out;
}

let F;   // кадр: { стиль: { lines, tris } }
const line = (style, pts) => {
    for (const r of runs(pts)) (STYLE[style].dash ? dashed(r, STYLE[style].dash) : [r]).forEach(l => F[style].lines.push(l));
};
const fan = (style, center, ring) => {
    for (let i = 1; i < ring.length; i++) F[style].tris.push(center[0], center[1], ring[i - 1][0], ring[i - 1][1], ring[i][0], ring[i][1]);
};
const disc = (style, c, r) => fan(style, c, Array.from({ length: 25 }, (_, i) => [c[0] + r * Math.cos(rad(i * 15)), c[1] + r * Math.sin(rad(i * 15))]));

function tessellate(lines, w) {
    const out = [], h = w / 2;
    for (const ln of lines) for (let i = 1; i < ln.length; i++) {
        const [ax, ay] = ln[i - 1], [bx, by] = ln[i], L = Math.hypot(bx - ax, by - ay); if (L < 1e-9) continue;
        const dx = (bx - ax) / L, dy = (by - ay) / L, nx = -dy * h, ny = dx * h;
        const x0 = ax - dx * h, y0 = ay - dy * h, x1 = bx + dx * h, y1 = by + dy * h;
        out.push(x0 + nx, y0 + ny, x0 - nx, y0 - ny, x1 + nx, y1 + ny, x1 + nx, y1 + ny, x0 - nx, y0 - ny, x1 - nx, y1 - ny);
    }
    return out;
}

/* ---------------------------------------------- 5. побудова кадру */
const niceStep = px => { const raw = px / S.view.sc, p = 10 ** Math.floor(Math.log10(raw)); return [1, 2, 5, 10].map(m => m * p).find(s => s >= raw) || p * 10; };

function buildFrame() {
    F = {}; for (const k of ORDER) F[k] = { lines: [], tris: [] };
    const T = eucl(), P = w => toPx(Mx.apply(T, w)), labels = [];
    const [x0, y0] = toWorld(0, 0), [x1, y1] = toWorld(W, H);

    // сітка та осі
    const st = niceStep(46);
    for (let i = Math.ceil(x0 / st); i <= x1 / st; i++) { const X = toPx([i * st, 0])[0]; line(i % 5 ? 'grid' : 'major', [[X, 0], [X, H]]); }
    for (let j = Math.ceil(y0 / st); j <= y1 / st; j++) { const Y = toPx([0, j * st])[1]; line(j % 5 ? 'grid' : 'major', [[0, Y], [W, Y]]); }
    const o = toPx([0, 0]);
    line('axisX', [[0, o[1]], [W, o[1]]]); line('axisY', [[o[0], 0], [o[0], H]]);
    for (let i = Math.ceil(x0 / st); i <= x1 / st; i++) if (i) labels.push([+(i * st).toFixed(4) + '', toPx([i * st, 0])[0] - 8, o[1] - 15, 'num']);
    for (let j = Math.ceil(y0 / st); j <= y1 / st; j++) if (j) labels.push([+(j * st).toFixed(4) + '', o[0] + 7, toPx([0, j * st])[1] - 4, 'num']);
    labels.push(['x', W - 336, o[1] + 8, 'axis'], ['y', o[0] + 9, H - 20, 'axis'], ['O', o[0] - 14, o[1] - 15, 'axis']);

    // крива: цикл за φ із кроком, щільнішим на кінцях гілок (де r → ∞)
    S.samples = [];
    const N = 2600;
    BRANCHES.forEach(([lo, hi], bi) => {
        const pts = [];
        for (let i = 0; i <= N; i++) {
            const phi = lo + 1e-5 + (hi - lo - 2e-5) * (1 - Math.cos(PI * i / N)) / 2;
            const px = P(K.pt(phi, S.a, S.b));
            pts.push(px);
            if (i % 3 === 0) S.samples.push({ phi, px, br: bi });
        }
        line('curve', pts);
    });
    // заливка петлі
    const lr = loopRange(S.a, S.b);
    if (S.show.loop && lr) {
        const ring = Array.from({ length: 241 }, (_, i) => P(K.pt(lr[0] + (lr[1] - lr[0]) * i / 240, S.a, S.b)));
        const cx = ring.reduce((s, p) => s + p[0], 0) / ring.length, cy = ring.reduce((s, p) => s + p[1], 0) / ring.length;
        fan('loop', [cx, cy], ring);
    }
    // асимптота x = a
    if (S.show.asy) {
        line('asy', [P([S.a, -1e4]), P([S.a, 1e4])]);
        const q = P([S.a, 0]); labels.push(['x = a', q[0] + 6, H - 40, 'asy']);
    }
    // точка M, дотична, нормаль, коло кривини
    const m = pointInfo();
    const half = 230;                                   // px в кожен бік — симетрична дотична
    const along = (dir, style) => {
        const c = P(m.w), d = Mx.apply(Mx.rotate(S.rot.a, 0, 0), dir);
        line(style, [[c[0] - d[0] * half, c[1] - d[1] * half], [c[0] + d[0] * half, c[1] + d[1] * half]]);
    };
    if (S.show.tan) along(m.tau, 'tan');
    if (S.show.nor) along([-m.tau[1], m.tau[0]], 'nor');
    if (S.show.circ && m.R * S.view.sc < 5000) {
        const s = Math.sign(m.k) || 1, c = P([m.w[0] - s * m.tau[1] * m.R, m.w[1] + s * m.tau[0] * m.R]);
        line('circ', Array.from({ length: 181 }, (_, i) => [c[0] + m.R * S.view.sc * Math.cos(rad(i * 2)), c[1] + m.R * S.view.sc * Math.sin(rad(i * 2))]));
        disc('dot', c, 2.5);
    }
    const mp = P(m.w); S.mPx = mp;
    disc('dot', mp, 6); labels.push(['M', mp[0] + 10, mp[1] + 9, 'mark']);
    if (S.rot.a) { const q = toPx([S.rot.x, S.rot.y]); line('pivot', [[q[0] - 10, q[1]], [q[0] + 10, q[1]]]); line('pivot', [[q[0], q[1] - 10], [q[0], q[1] + 10]]); }
    return labels;
}

function pointInfo() {
    const phi = rad(S.phi), a = S.a, b = S.b, w = K.pt(phi, a, b), d = K.d1(phi, a, b), v = Math.hypot(d[0], d[1]), k = K.curv(phi, a, b);
    return { phi, w, tau: [d[0] / v, d[1] / v], k, R: 1 / Math.abs(k) };
}

function paint() {
    const labels = buildFrame();
    g2.setTransform(dpr, 0, 0, dpr, 0, 0); g2.clearRect(0, 0, W, H);
    if (gl) {
        gl.viewport(0, 0, glc.width, glc.height);
        gl.clearColor(0.047, 0.165, 0.302, 1); gl.clear(gl.COLOR_BUFFER_BIT);
        gl.uniform2f(loc.res, W, H); gl.bindBuffer(gl.ARRAY_BUFFER, buf);
        gl.enableVertexAttribArray(loc.p); gl.vertexAttribPointer(loc.p, 2, gl.FLOAT, false, 0, 0);
        for (const k of ORDER) {
            const tris = FILLS.has(k) ? F[k].tris : tessellate(F[k].lines, STYLE[k].w);
            if (!tris.length) continue;
            gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(tris), gl.DYNAMIC_DRAW);
            gl.uniform4fv(loc.col, STYLE[k].c); gl.drawArrays(gl.TRIANGLES, 0, tris.length / 2);
        }
    } else {
        g2.fillStyle = '#0c2a4d'; g2.fillRect(0, 0, W, H);
        for (const k of ORDER) {
            const s = STYLE[k], col = `rgba(${s.c.slice(0, 3).map(v => v * 255 | 0)},${s.c[3]})`;
            g2.strokeStyle = g2.fillStyle = col; g2.lineWidth = s.w; g2.lineJoin = 'round';
            if (FILLS.has(k)) { const t = F[k].tris; for (let i = 0; i < t.length; i += 6) { g2.beginPath(); g2.moveTo(t[i], H - t[i + 1]); g2.lineTo(t[i + 2], H - t[i + 3]); g2.lineTo(t[i + 4], H - t[i + 5]); g2.fill(); } }
            else for (const ln of F[k].lines) { g2.beginPath(); ln.forEach((p, i) => i ? g2.lineTo(p[0], H - p[1]) : g2.moveTo(p[0], H - p[1])); g2.stroke(); }
        }
    }
    for (const [t, x, y, kind] of labels) {
        g2.font = kind === 'mark' ? 'italic 700 15px "Times New Roman", serif' : kind === 'axis' ? 'italic 14px "Times New Roman", serif' : kind === 'asy' ? 'italic 13px Cambria, serif' : '10.5px Consolas, monospace';
        g2.fillStyle = kind === 'mark' ? '#ffd166' : kind === 'asy' ? '#ff8f8f' : kind === 'axis' ? '#d1ebff' : '#6f9bc9';
        g2.fillText(t, x, H - y);
    }
    report();
}

/* ------------------------------------------------------- 6. звіт */
const f4 = v => Number.isFinite(v) ? (Math.round(v * 1e4) / 1e4).toString() : '∞';
const f5 = v => Number.isFinite(v) ? v.toFixed(5) : '∞';
const wrap360 = d => ((d + 90) % 360 + 360) % 360 - 90;
function report() {
    const a = S.a, b = S.b, m = pointInfo(), T = eucl(), Mw = Mx.apply(T, m.w);
    const [p1, p2] = S.calc.map(rad), L = arcLength(p1, p2, a, b), lr = loopRange(a, b);
    const infl = inflections(a, b);
    const lines = [
        `M: φ = ${f4(S.phi)}°   (${f4(m.w[0])}; ${f4(m.w[1])})`,
        `після перетворень: (${f4(Mw[0])}; ${f4(Mw[1])})`,
        `k = ${f5(m.k)},   R = ${f5(m.R)}`,
        '',
        `L[${S.calc[0]}°; ${S.calc[1]}°] = ${Number.isFinite(L) ? f5(L) : '∞ (полюс cos φ = 0)'}`,
        lr ? `Петля: φ ∈ [${f4(deg(lr[0]))}°; ${f4(deg(lr[1]))}°]` : 'Петлі немає (b ≤ a)',
        lr ? `  довжина = ${f5(arcLength(lr[0], lr[1], a, b))}` : '',
        lr ? `  площа (числ.) = ${f5(loopAreaNumeric(a, b))}` : '',
        lr ? `  площа (форм.) = ${f5(loopAreaExact(a, b))}` : '',
        `Площа між віткою та асимптотою = ∞`,
        '',
        `Асимптота: x = ${f4(a)}`,
        `Перегини (${infl.length}):`,
        ...infl.map(p => { const w = K.pt(p, a, b); return `  φ = ${f4(wrap360(deg(p)))}°  (${f4(w[0])}; ${f4(w[1])})`; })
    ].filter((s, i, arr) => s !== '' || (arr[i - 1] !== '' && i > 0));
    $('calc').textContent = lines.join('\n');
    $('aO').textContent = f4(S.a); $('bO').textContent = f4(S.b); $('phiO').textContent = f4(S.phi) + '°'; $('raO').textContent = S.rot.a + '°';
}

/* ---------------------------------------------- 7. інтерфейс */
function syncInputs() {
    $('a').value = S.a; $('b').value = S.b; $('phi').value = S.phi; $('ra').value = S.rot.a;
    $('shx').value = S.shift[0]; $('shy').value = S.shift[1]; $('rx').value = S.rot.x; $('ry').value = S.rot.y;
    $('c1').value = S.calc[0]; $('c2').value = S.calc[1];
}
function animTick(now) {
    if (!S.anim.on) return;
    const t = (now - S.anim.t0) / 1000 * S.anim.speed * 0.25;
    S.a = Math.round((2.3 + 1.2 * Math.sin(t)) * 100) / 100;
    S.b = Math.round((3.4 + 1.9 * Math.sin(t * 0.61 + 1)) * 100) / 100;
    syncInputs(); paint(); requestAnimationFrame(animTick);
}
function setAnim(on) {
    S.anim.on = on; $('btnAnim').textContent = on ? '■ Stop' : '▶ Animate';
    if (on) { S.anim.t0 = performance.now(); requestAnimationFrame(animTick); }
}
const nearestIdx = () => { let bi = 0, bd = 1e18; S.samples.forEach((s, i) => { const d = Math.abs(deg(s.phi) - S.phi); if (d < bd) { bd = d; bi = i; } }); return bi; };
const sorted = () => S.samples.map(s => deg(s.phi)).sort((u, v) => u - v);

function dragPoint(mx, my) {
    // найближча вибірка в екранних координатах + уточнення тернарним пошуком у межах гілки
    let best = null, bd = 1e18;
    for (const s of S.samples) { const d = (s.px[0] - mx) ** 2 + (s.px[1] - my) ** 2; if (d < bd) { bd = d; best = s; } }
    if (!best) return;
    const [lo0, hi0] = BRANCHES[best.br], T = eucl();
    const dist = p => { const q = toPx(Mx.apply(T, K.pt(p, S.a, S.b))); return (q[0] - mx) ** 2 + (q[1] - my) ** 2; };
    let lo = Math.max(lo0 + 1e-3, best.phi - 0.08), hi = Math.min(hi0 - 1e-3, best.phi + 0.08);
    for (let i = 0; i < 40; i++) { const m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3; dist(m1) < dist(m2) ? hi = m2 : lo = m1; }
    S.phi = Math.round(deg((lo + hi) / 2) * 100) / 100; syncInputs(); paint();
}

function bind() {
    const slider = (id, fn) => $(id).oninput = () => { fn(parseFloat($(id).value)); paint(); };
    slider('a', v => S.a = v); slider('b', v => S.b = v); slider('phi', v => S.phi = v);
    slider('ra', v => S.rot.a = v);
    $('speed').oninput = () => { S.anim.speed = +$('speed').value; $('spO').textContent = S.anim.speed; };
    $('btnDefault').onclick = () => { S.a = 2; S.b = 3; S.phi = 40; S.view.sc = 42; syncInputs(); paint(); };
    $('btnAnim').onclick = () => setAnim(!S.anim.on);
    for (const [id, k] of [['vTan', 'tan'], ['vNor', 'nor'], ['vAsy', 'asy'], ['vCirc', 'circ'], ['vLoop', 'loop']]) $(id).onchange = e => { S.show[k] = e.target.checked; paint(); };
    $('pPrev').onclick = () => { const l = sorted(), c = S.phi; S.phi = l.filter(v => v < c - 0.01).pop() ?? l[0]; syncInputs(); paint(); };
    $('pNext').onclick = () => { const l = sorted(), c = S.phi; S.phi = l.find(v => v > c + 0.01) ?? l[l.length - 1]; syncInputs(); paint(); };
    for (const id of ['shx', 'shy', 'rx', 'ry']) $(id).oninput = () => { S.shift = [num('shx'), num('shy')]; S.rot.x = num('rx'); S.rot.y = num('ry'); paint(); };
    $('btnReset').onclick = () => { S.shift = [0, 0]; S.rot = { x: 0, y: 0, a: 0 }; syncInputs(); paint(); };
    $('c1').oninput = $('c2').oninput = () => { S.calc = [num('c1'), num('c2')]; paint(); };
    $('presetLoop').onclick = () => { const lr = loopRange(S.a, S.b); S.calc = lr ? [+deg(lr[0]).toFixed(3), +deg(lr[1]).toFixed(3)] : [150, 210]; syncInputs(); paint(); };
    $('presetAll').onclick = () => { S.calc = [-60, 60]; syncInputs(); paint(); };
    document.querySelectorAll('[data-fold]').forEach(h => h.onclick = () => h.parentElement.classList.toggle('folded'));

    // миша: панорамування, перетягування M, масштаб до курсору
    let mode = null, start = null;
    const pos = e => [e.clientX, H - e.clientY];
    glc.addEventListener('pointerdown', e => {
        const [x, y] = pos(e);
        if (e.shiftKey) { const w = toWorld(x, y); S.rot.x = Math.round(w[0] * 4) / 4; S.rot.y = Math.round(w[1] * 4) / 4; if (!S.rot.a) S.rot.a = 0.0001; syncInputs(); paint(); return; }
        mode = S.mPx && Math.hypot(S.mPx[0] - x, S.mPx[1] - y) < 16 ? 'M' : 'pan';
        start = { x, y, ox: S.view.ox, oy: S.view.oy };
        glc.setPointerCapture(e.pointerId); glc.classList.add('drag');
    });
    glc.addEventListener('pointermove', e => {
        if (!mode) return; const [x, y] = pos(e);
        if (mode === 'M') dragPoint(x, y);
        else { S.view.ox = start.ox + x - start.x; S.view.oy = start.oy + y - start.y; paint(); }
    });
    glc.addEventListener('pointerup', () => { mode = null; glc.classList.remove('drag'); });
    glc.addEventListener('wheel', e => {
        e.preventDefault();
        const [x, y] = pos(e), k = e.deltaY < 0 ? 1.12 : 1 / 1.12, sc = Math.min(160, Math.max(8, S.view.sc * k)), r = sc / S.view.sc;
        S.view.ox = x - (x - S.view.ox) * r; S.view.oy = y - (y - S.view.oy) * r; S.view.sc = sc; paint();
    }, { passive: false });
    window.addEventListener('keydown', e => { if (e.target.tagName !== 'INPUT' && (e.key === 'a' || e.key === 'A')) setAnim(!S.anim.on); });
    window.addEventListener('resize', layout);
}

let firstLayout = true;
function layout() {
    dpr = window.devicePixelRatio || 1; W = innerWidth; H = innerHeight;
    for (const c of [glc, ov]) { c.width = Math.round(W * dpr); c.height = Math.round(H * dpr); }
    if (firstLayout) { S.view.ox = W * 0.40; S.view.oy = H * 0.46; firstLayout = false; }
    paint();
}

function applyQuery() {
    const q = new URLSearchParams(location.search).get('s'); if (!q) return;
    try {
        const s = JSON.parse(q);
        for (const k of ['a', 'b', 'phi']) if (s[k] !== undefined) S[k] = s[k];
        if (s.shift) S.shift = s.shift;
        if (s.rot) Object.assign(S.rot, s.rot);
        if (s.sc) S.view.sc = s.sc;
        if (s.circ) { S.show.circ = true; $('vCirc').checked = true; }
        if (s.calc) S.calc = s.calc;
        syncInputs();
        if (s.anim) setAnim(true);
    } catch (e) { /* ignore */ }
}

initGL();
$('renderer').textContent = gl ? 'WebGL' : 'Canvas 2D (резервний режим)';
bind(); applyQuery(); layout();
