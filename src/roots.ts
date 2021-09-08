const GOLDEN_RATIO = (1 + Math.sqrt(5)) / 2;

const getMachineEpsilon = (): number => {
  let e = 1.0;
  while (1.0 + e !== 1.0) {
    e *= 0.5;
  }
  return e;
};

const MACHINE_EPSILON: number = getMachineEpsilon();

type NumToNumFn = (x: number) => number;

const sign = (x: number): number => {
  if (x) {
    if (x < 0) {
      return -1;
    } else {
      return 1;
    }
  } else {
    if (x === x) {
      return 0;
    } else {
      return NaN;
    }
  }
};

// Finds the root of f(x) near x0 given df(x) = f'(x)
const newtonsMethod = (x0: number, f: NumToNumFn, df: NumToNumFn): number => {
  let x: number;
  while (true) {
    x = x0 - f(x0) / df(x0);
    if (isNaN(x) || Math.abs(x - x0) < 1e-6) {
      return x;
    }
    x0 = x;
  }
};

// Finds a root of f(x) between a and b
const brentsMethod = (a: number, b: number, relativeAccuracy: number, f: NumToNumFn, fa?: number, fb?: number) => {
  let c: number;
  let d: number;
  let e: number;
  let fc: number;
  let m: number;
  let p: number;
  let q: number;
  let r: number;
  let s: number;
  let tol: number;

  if (fa == null) {
    fa = f(a);
  }
  if (fb == null) {
    fb = f(b);
  }

  c = a;
  fc = fa;
  d = b - a;
  e = d;
  relativeAccuracy += 0.5 * MACHINE_EPSILON;

  if (isNaN(fa) || isNaN(fb)) {
    return NaN;
  }
  if (sign(fa) === sign(fb)) {
    // Can't find a root if the signs of fa and fb are equal
    return NaN;
  }

  for (let i = 0; i < 100; i++) {
    if (Math.abs(fc) < Math.abs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol = relativeAccuracy * Math.abs(b);
    m = 0.5 * (c - b);
    if (fb === 0 || Math.abs(m) <= tol) {
      return b;
    }
    if (Math.abs(e) < tol || Math.abs(fa) <= Math.abs(fb)) {
      d = e = m;
    } else {
      s = fb / fa;
      if (a === c) {
        // Use a linear interpolation step
        p = 2 * m * s;
        q = 1 - s;
      } else {
        // Use a parabolic interpolation step
        q = fa / fc;
        r = fb / fc;
        p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
        q = (q - 1) * (r - 1) * (s - 1);
      }
      if (p > 0) {
        q = -q;
      } else {
        p = -p;
      }
      // Validate interpolation
      if (2 * p < Math.min(3 * m * q - Math.abs(tol * q), Math.abs(e * q))) {
        e = d;
        d = p / q;
      } else {
        // Fall back to bisection
        d = e = m;
      }
    }
    a = b;
    fa = fb;
    if (Math.abs(d) > tol) {
      b += d;
    } else {
      b += m > 0 ? tol : -tol;
    }
    fb = f(b);
    if (isNaN(fb)) {
      return NaN;
    }
    if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) {
      c = a;
      fc = fa;
      d = e = b - a;
    }
  }
  throw Error("Brent's method failed to converge after 100 iterations");
};

// Finds the minimum of f(x) between x1 and x2. Returns x.
// See: http://en.wikipedia.org/wiki/Golden_section_search
const goldenSectionSearch = (x1: number, x2: number, epsilon: number, f: NumToNumFn): number => {
  let k: number;
  let x: number;
  let x3: number;
  let y: number;
  let y2: number;

  k = 2 - GOLDEN_RATIO;
  x3 = x2;
  x2 = x1 + k * (x3 - x1);
  y2 = f(x2);
  while (true) {
    if (x3 - x2 > x2 - x1) {
      x = x2 + k * (x3 - x2);
    } else {
      x = x2 - k * (x2 - x1);
    }
    if (x3 - x1 < epsilon * (x2 + x)) {
      return (x3 + x1) / 2;
    }
    y = f(x);
    if (y < y2) {
      if (x3 - x2 > x2 - x1) {
        x1 = x2;
      } else {
        x3 = x2;
      }
      x2 = x;
      y2 = y;
    } else {
      if (x3 - x2 > x2 - x1) {
        x3 = x;
      } else {
        x1 = x;
      }
    }
  }
};

export { MACHINE_EPSILON, newtonsMethod, brentsMethod, goldenSectionSearch };
