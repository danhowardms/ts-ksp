import { MACHINE_EPSILON, brentsMethod } from "./roots";
import { addVV, dotVV, mulVS, normV, subVV, angleBetween, Vector3 } from "./vector3";

const TWO_PI = 2 * Math.PI;
const HALF_PI = 0.5 * Math.PI;

const acot = (x: number): number => {
  return HALF_PI - Math.atan(x);
};

const acoth = (x: number): number => {
  return 0.5 * Math.log((x + 1) / (x - 1));
};

const relativeError = (a: number, b: number): number => {
  return Math.abs(1.0 - a / b);
};

type LambertSolution = {
  ejectionVelocity: Vector3;
  insertionVelocity: Vector3;
  transferAngle: number;
};

/**
 * Find solutions to Lambert's problem:
 * Find the orbit that passes through two given points at two given times
 *
 * @param mu Rate of acceleration towards the central body (m/s2 I think)
 * @param pos1 Position vector at the start of the transfer
 * @param pos2 Position vector at the end of the transfer
 * @param dt Time to move from pos1 to pos2
 * @param maxRevs Max number of revolutions around an elliptical transfer orbit to consider
 * @param prograde ??? I'm not sure what this parameter means exactly, but it's always omitted where it was used in
 *        AlexMoon's transfer calculator, which means it would always take the default value, 1. Maybe -1 would give
 *        us a retrograde transfer?
 */
const solveLambert = (mu: number, pos1: Vector3, pos2: Vector3, dt: number, maxRevs?: number, prograde?: number): LambertSolution[] => {
  // Based on Sun, F.T. "On the Minimum Time Trajectory and Multiple Solutions of Lambert's Problem"
  // AAS/AIAA Astrodynamics Conference, Provincetown, Massachusetts, AAS 79-164, June 25-27, 1979

  if (maxRevs == null) {
    maxRevs = 0;
  }
  if (prograde == null) {
    prograde = 1;
  }

  let N = 0;
  const r1 = normV(pos1);
  const r2 = normV(pos2);

  // Intermediate terms
  const deltaPos = subVV(pos2, pos1);
  const c = normV(deltaPos);
  const m = r1 + r2 + c;
  const n = r1 + r2 - c;

  // Assume we want a prograde orbit counter-clockwise around the +z axis
  let transferAngle = angleBetween(pos1, pos2);
  if ((pos1[0] * pos2[1] - pos1[1] * pos2[0]) * prograde < 0) {
    transferAngle = TWO_PI - transferAngle;
  }

  let angleParameter = Math.sqrt(n / m);
  if (transferAngle > Math.PI) {
    angleParameter = -angleParameter;
  }

  const normalizedTime = 4 * dt * Math.sqrt(mu / (m * m * m));
  const parabolicNormalizedTime = (2 / 3) * (1 - angleParameter * angleParameter * angleParameter);

  // Pre-calculate terms for efficiency
  const sqrtMu = Math.sqrt(mu);
  const invSqrtM = 1 / Math.sqrt(m);
  const invSqrtN = 1 / Math.sqrt(n);

  const solutions: LambertSolution[] = [];
  const pushSolution = (_x: number, _y: number, _N: number): void => {
    const vc = sqrtMu * (_y * invSqrtN + _x * invSqrtM);
    const vr = sqrtMu * (_y * invSqrtN - _x * invSqrtM);
    const ec = mulVS(deltaPos, vc / c);
    const v1 = addVV(ec, mulVS(pos1, vr / r1));
    const v2 = subVV(ec, mulVS(pos2, vr / r2));
    solutions.push({
      ejectionVelocity: v1,
      insertionVelocity: v2,
      transferAngle: _N * TWO_PI + transferAngle,
    });
  };

  // y = +/- sqrt(1 - sigma^2 * (1 - x^2))
  const fy = (_x: number): number => {
    const _y = Math.sqrt(1 - angleParameter * angleParameter * (1 - _x * _x));
    if (angleParameter < 0) {
      return -_y;
    } else {
      return _y;
    }
  };

  // Returns the difference our desired normalizedTime and the normalized
  // time for a path parameter of x (given our angleParameter)
  // Defined over the domain of (-1, infinity)
  const ftau = (_x: number): number => {
    if (_x === 1.0) {
      // Parabolic
      return parabolicNormalizedTime - normalizedTime;
    } else {
      const _y = fy(_x);
      if (_x > 1) {
        // Hyperbolic
        const g = Math.sqrt(_x * _x - 1);
        const h = Math.sqrt(_y * _y - 1);
        return (-acoth(_x / g) + acoth(_y / h) + _x * g - _y * h) / (g * g * g) - normalizedTime;
      } else {
        // Elliptical: -1 < x < 1
        const g = Math.sqrt(1 - _x * _x);
        const h = Math.sqrt(1 - _y * _y);
        return (acot(_x / g) - Math.atan(h / _y) - _x * g + _y * h + N * Math.PI) / (g * g * g) - normalizedTime;
      }
    }
  };

  let x: number;
  let y: number;
  let x1: number;
  let x2: number;
  let xMT: number;
  let minimumNormalizedTime: number;

  // Partition the solution space
  if (relativeError(normalizedTime, parabolicNormalizedTime) < 1e-6) {
    // Unique parabolic solution
    x = 1.0;
    y = angleParameter < 0 ? -1 : 1;
    pushSolution(x, y, 0);
  } else if (normalizedTime < parabolicNormalizedTime) {
    // Unique hyperbolic solution
    x1 = 1.0;
    x2 = 2.0;
    while (!(ftau(x2) < 0.0)) {
      x1 = x2;
      x2 *= 2.0;
      if (! isFinite(x2)) {
        throw Error("Cannot solve Lambert for hyperbolic - x got too big");
      }
    }
    x = brentsMethod(x1, x2, 1e-4, ftau);
    pushSolution(x, fy(x), N);
  } else {
    // Potentially multiple elliptical solutions
    maxRevs = Math.min(maxRevs, Math.floor(normalizedTime / Math.PI));
    let minimumEnergyNormalizedTime = Math.acos(angleParameter) + angleParameter * Math.sqrt(1 - angleParameter * angleParameter);

    for (N = 0; N <= maxRevs; N++) {
      if (N > 0 && N === maxRevs) {
        // Check the number of solutions for the last revolution
        const phix = (_x: number) => {
          const g = Math.sqrt(1 - _x * _x);
          return acot(_x / g) - ((2 + _x * _x) * g) / (3 * _x);
        };

        const phiy = (_y: number) => {
          const h = Math.sqrt(1 - _y * _y);
          return Math.atan(h / _y) - ((2 + _y * _y) * h) / (3 * _y);
        };

        // Find the minimum (normalized) time an N revolution trajectory will take
        if (angleParameter === 1) {
          xMT = 0;
          minimumNormalizedTime = minimumEnergyNormalizedTime;
        } else if (angleParameter === 0) {
          xMT = brentsMethod(0, 1, 1e-6, (_x: number): number => {
            return phix(_x) + N * Math.PI;
          });
          minimumNormalizedTime = 2 / (3 * xMT);
        } else {
          xMT = brentsMethod(0, 1, 1e-6, (_x: number): number => {
            return phix(_x) - phiy(fy(_x)) + N * Math.PI;
          });
          minimumNormalizedTime = (2 / 3) * (1 / xMT - (angleParameter * angleParameter * angleParameter) / Math.abs(fy(xMT)));
        }

        if (relativeError(normalizedTime, minimumNormalizedTime) < 1e-6) {
          // One solution for N revolutions and we're done
          pushSolution(xMT, fy(xMT), (N + 1) * TWO_PI - transferAngle);
          break;
        } else if (normalizedTime < minimumNormalizedTime) {
          // No solutions for N revolutions; we're done
          break;
        } else if (normalizedTime < minimumEnergyNormalizedTime) {
          // Two low path solutions
          x = brentsMethod(0, xMT, 1e-4, ftau);
          if (!isNaN(x)) {
            pushSolution(x, fy(x), N);
          }
          x = brentsMethod(xMT, 1.0 - MACHINE_EPSILON, 1e-4, ftau);
          if (!isNaN(x)) {
            pushSolution(x, fy(x), N);
          }
          break;
        }
      }

      if (relativeError(normalizedTime, minimumEnergyNormalizedTime) < 1e-6) {
        // The minimum energy elliptical solution
        pushSolution(0, fy(0), N);
        if (N > 0) {
          x = brentsMethod(1e-6, 1.0 - MACHINE_EPSILON, 1e-4, ftau);
          if (!isNaN(x)) {
            pushSolution(x, fy(x), N);
          }
        }
      } else {
        if (N > 0 || normalizedTime > minimumEnergyNormalizedTime) {
          // High path solution
          x = brentsMethod(-1.0 + MACHINE_EPSILON, 0, 1e-4, ftau);
          if (!isNaN(x)) {
            pushSolution(x, fy(x), N);
          }
        }
        if (N > 0 || normalizedTime < minimumEnergyNormalizedTime) {
          // Low path solution
          x = brentsMethod(0, 1.0 - MACHINE_EPSILON, 1e-4, ftau);
          if (!isNaN(x)) {
            pushSolution(x, fy(x), N);
          }
        }
      }
      minimumEnergyNormalizedTime += Math.PI;
    }
  }
  return solutions;
};

export { LambertSolution, solveLambert };
