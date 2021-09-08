import { Vector3, addVV, subVV, mulVS, divVS, normSquaredV, normV, normalizeV, dotVV, crossVV } from "./vector3";
import { Quaternion, conjugateQ, concatQQ, quaternionFromAngleAndAxis, quaternionFromStartAndEndVectors, rotate } from "./quaternion";
import { newtonsMethod, goldenSectionSearch } from "./roots";
import { solveLambert } from "./lambert";
import { TransferOptions } from "./transfer-options";
import { CelestialBody } from "./celestial-body";
import { AngleDegrees } from "./utility-types";
import { OrbitingCelestialBody } from "./orbiting-celestial-body";

const TWO_PI: number = 2 * Math.PI;
const HALF_PI: number = 0.5 * Math.PI;

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

const sinh = (angle: number): number => {
  const p = Math.exp(angle);
  return (p - 1 / p) * 0.5;
};

const cosh = (angle: number): number => {
  const p = Math.exp(angle);
  return (p + 1 / p) * 0.5;
};

const acosh = (n: number): number => {
  return Math.log(n + Math.sqrt(n * n - 1));
};

// Find the projection of the point p onto the plane (through the origin) defined by the unit normal vector n
const projectToPlane = (p: Vector3, n: Vector3): Vector3 => {
  return subVV(p, mulVS(n, dotVV(p, n)));
};

// @todo I THINK what is happening here is that from and to are both projected into the plane defined by normal
// Then we find out the angle between their projections
const angleInPlane = (from: Vector3, to: Vector3, normal: Vector3): number => {
  from = normalizeV(projectToPlane(from, normal));
  to = normalizeV(projectToPlane(to, normal));
  const rot = quaternionFromStartAndEndVectors(normal, [0, 0, 1]);
  from = rotate(rot, from);
  to = rotate(rot, to);
  const result = Math.atan2(from[1], from[0]) - Math.atan2(to[1], to[0]);
  if (result < 0) {
    return result + TWO_PI;
  } else {
    return result;
  }
};

const degreesToRadians = (degrees: AngleDegrees): number => {
  return ((degrees * Math.PI) / 180) as number;
};

class Orbit {
  referenceBody: CelestialBody;
  semiMajorAxis: number;
  eccentricity: number;
  inclination: number = 0;
  longitudeOfAscendingNode: number = 0;
  argumentOfPeriapsis: number = 0;
  meanAnomalyAtEpoch: number = 0;
  timeOfPeriapsisPassage: number = 0;

  constructor(referenceBody1: CelestialBody, semiMajorAxis1: number, eccentricity1: number, inclinationDegs?: AngleDegrees, longitudeOfAscendingNodeDegs?: AngleDegrees, argumentOfPeriapsisDegs?: AngleDegrees, meanAnomalyAtEpoch1?: number, timeOfPeriapsisPassage1?: number) {
    this.referenceBody = referenceBody1;
    this.semiMajorAxis = semiMajorAxis1;
    this.eccentricity = eccentricity1;
    if (inclinationDegs != null) {
      this.inclination = degreesToRadians(inclinationDegs);
    }
    if (longitudeOfAscendingNodeDegs != null) {
      this.longitudeOfAscendingNode = degreesToRadians(longitudeOfAscendingNodeDegs);
    }
    if (argumentOfPeriapsisDegs != null) {
      this.argumentOfPeriapsis = degreesToRadians(argumentOfPeriapsisDegs);
    }
    if (meanAnomalyAtEpoch1 != null) {
      this.meanAnomalyAtEpoch = meanAnomalyAtEpoch1;
    }
    if (timeOfPeriapsisPassage1 != null) {
      this.timeOfPeriapsisPassage = timeOfPeriapsisPassage1;
    }
  }

  isHyperbolic(): boolean {
    return this.eccentricity > 1;
  }

  apoapsis(): number {
    return this.semiMajorAxis * (1 + this.eccentricity);
  }

  periapsis(): number {
    return this.semiMajorAxis * (1 - this.eccentricity);
  }

  apoapsisAltitude(): number {
    return this.apoapsis() - this.referenceBody.radius;
  }

  periapsisAltitude(): number {
    return this.periapsis() - this.referenceBody.radius;
  }

  semiMinorAxis(): number {
    const e = this.eccentricity;
    return this.semiMajorAxis * Math.sqrt(1 - e * e);
  }

  semiLatusRectum(): number {
    const e = this.eccentricity;
    return this.semiMajorAxis * (1 - e * e);
  }

  meanMotion(): number {
    const a = Math.abs(this.semiMajorAxis);
    return Math.sqrt(this.referenceBody.gravitationalParameter / (a * a * a));
  }

  period(): number {
    if (this.isHyperbolic()) {
      return Infinity;
    } else {
      return TWO_PI / this.meanMotion();
    }
  }

  rotationToReferenceFrame(): Quaternion {
    const axisOfInclination: Vector3 = [Math.cos(-this.argumentOfPeriapsis), Math.sin(-this.argumentOfPeriapsis), 0];
    return concatQQ(quaternionFromAngleAndAxis(this.longitudeOfAscendingNode + this.argumentOfPeriapsis, [0, 0, 1]), quaternionFromAngleAndAxis(this.inclination, axisOfInclination));
  }

  normalVector(): Vector3 {
    return rotate(this.rotationToReferenceFrame(), [0, 0, 1]);
  }

  // @todo is orbit a different orbit?
  phaseAngle(orbit: Orbit, t: number): number {
    const n: Vector3 = this.normalVector();
    const p1: Vector3 = this.positionAtTrueAnomaly(this.trueAnomalyAt(t));
    let p2: Vector3 = orbit.positionAtTrueAnomaly(orbit.trueAnomalyAt(t));
    p2 = subVV(p2, mulVS(n, dotVV(p2, n)));
    const r1: number = normV(p1);
    const r2: number = normV(p2);
    let phaseAngle: number = Math.acos(dotVV(p1, p2) / (r1 * r2));
    if (dotVV(crossVV(p1, p2), n) < 0) {
      phaseAngle = TWO_PI - phaseAngle;
    }
    if (orbit.semiMajorAxis < this.semiMajorAxis) {
      phaseAngle = phaseAngle - TWO_PI;
    }
    return phaseAngle;
  }

  meanAnomalyAt(t: number): number {
    if (this.isHyperbolic()) {
      return (t - this.timeOfPeriapsisPassage) * this.meanMotion();
    } else {
      if (this.timeOfPeriapsisPassage != null) {
        const M = ((t - this.timeOfPeriapsisPassage) % this.period()) * this.meanMotion();
        if (M < 0) {
          return M + TWO_PI;
        } else {
          return M;
        }
      } else {
        return (this.meanAnomalyAtEpoch + this.meanMotion() * (t % this.period())) % TWO_PI;
      }
    }
  }

  eccentricAnomalyAt(t: number): number {
    const e = this.eccentricity;
    const M = this.meanAnomalyAt(t);
    if (this.isHyperbolic()) {
      return newtonsMethod(M, (x) => {
        return M - e * sinh(x) + x;
      }, (x) => {
        return 1 - e * cosh(x);
      });
    } else {
      return newtonsMethod(M, (x) => {
        return M + e * Math.sin(x) - x;
      }, (x) => {
        return e * Math.cos(x) - 1;
      });
    }
  }

  trueAnomalyAt(t: number): number {
    const e = this.eccentricity;
    if (this.isHyperbolic()) {
      const H = this.eccentricAnomalyAt(t);
      const tA = Math.acos((e - cosh(H)) / (cosh(H) * e - 1));
      if (H < 0) {
        return -tA;
      } else {
        return tA;
      }
    } else {
      const E = this.eccentricAnomalyAt(t);
      const tA = 2 * Math.atan2(Math.sqrt(1 + e) * Math.sin(E / 2), Math.sqrt(1 - e) * Math.cos(E / 2));
      if (tA < 0) {
        return tA + TWO_PI;
      } else {
        return tA;
      }
    }
  }

  eccentricAnomalyAtTrueAnomaly(tA: number): number {
    const e = this.eccentricity;
    if (this.isHyperbolic()) {
      const cosTrueAnomaly = Math.cos(tA);
      const H = acosh((e + cosTrueAnomaly) / (1 + e * cosTrueAnomaly));
      if (tA < 0) {
        return -H;
      } else {
        return H;
      }
    } else {
      const E = 2 * Math.atan(Math.tan(tA / 2) / Math.sqrt((1 + e) / (1 - e)));
      if (E < 0) {
        return E + TWO_PI;
      } else {
        return E;
      }
    }
  }

  meanAnomalyAtTrueAnomaly(tA: number): number {
    const e = this.eccentricity;
    if (this.isHyperbolic()) {
      const H = this.eccentricAnomalyAtTrueAnomaly(tA);
      return e * sinh(H) - H;
    } else {
      const E = this.eccentricAnomalyAtTrueAnomaly(tA);
      return E - e * Math.sin(E);
    }
  }

  timeAtTrueAnomaly(tA: number, t0: number | null): number {
    if (t0 == null) {
      t0 = 0;
    }
    const M = this.meanAnomalyAtTrueAnomaly(tA);
    if (this.isHyperbolic()) {
      return this.timeOfPeriapsisPassage + M / this.meanMotion();
    } else {
      let t: number;
      const p = this.period();
      if (this.timeOfPeriapsisPassage != null) {
        t = this.timeOfPeriapsisPassage + p * Math.floor((t0 - this.timeOfPeriapsisPassage) / p) + M / this.meanMotion();
      } else {
        t = t0 - (t0 % p) + (M - this.meanAnomalyAtEpoch) / this.meanMotion();
      }
      if (t < t0) {
        return t + p;
      } else {
        return t;
      }
    }
  }

  radiusAtTrueAnomaly(tA: number): number {
    const e = this.eccentricity;
    return (this.semiMajorAxis * (1 - e * e)) / (1 + e * Math.cos(tA));
  }

  altitudeAtTrueAnomaly(tA: number): number {
    return this.radiusAtTrueAnomaly(tA) - this.referenceBody.radius;
  }

  speedAtTrueAnomaly(tA: number): number {
    return Math.sqrt(this.referenceBody.gravitationalParameter * (2 / this.radiusAtTrueAnomaly(tA) - 1 / this.semiMajorAxis));
  }

  positionAtTrueAnomaly(tA: number): Vector3 {
    const r = this.radiusAtTrueAnomaly(tA);
    return rotate(this.rotationToReferenceFrame(), [r * Math.cos(tA), r * Math.sin(tA), 0]);
  }

  velocityAtTrueAnomaly(tA: number): Vector3 {
    const mu = this.referenceBody.gravitationalParameter;
    const e = this.eccentricity;
    const h = Math.sqrt(mu * this.semiMajorAxis * (1 - e * e));
    const r = this.radiusAtTrueAnomaly(tA);
    const sin = Math.sin(tA);
    const cos = Math.cos(tA);
    const vr = (mu * e * sin) / h;
    const vtA = h / r;
    return rotate(this.rotationToReferenceFrame(), [vr * cos - vtA * sin, vr * sin + vtA * cos, 0]);
  }

  trueAnomalyAtPosition(p: Vector3): number {
    p = rotate(conjugateQ(this.rotationToReferenceFrame()), p);
    return Math.atan2(p[1], p[0]);
  }

  /* @todo
    static fromJSON(json): Orbit {
        const referenceBody = CelestialBody.fromJSON(json.referenceBody);
        const orbit = new Orbit(referenceBody, json.semiMajorAxis, json.eccentricity);
        orbit.inclination = json.inclination;
        orbit.longitudeOfAscendingNode = json.longitudeOfAscendingNode;
        orbit.argumentOfPeriapsis = json.argumentOfPeriapsis;
        orbit.meanAnomalyAtEpoch = json.meanAnomalyAtEpoch;
        orbit.timeOfPeriapsisPassage = json.timeOfPeriapsisPassage;
        return orbit;
    };
     */

  static fromApoapsisAndPeriapsis(referenceBody: CelestialBody, apoapsis: number, periapsis: number, inclinationDegs: AngleDegrees, longitudeOfAscendingNodeDegs: AngleDegrees, argumentOfPeriapsisDegs: AngleDegrees, meanAnomalyAtEpoch: number, timeOfPeriapsisPassage: number) {
    if (apoapsis < periapsis) {
      // swap them over
      const tmp: number = apoapsis;
      apoapsis = periapsis;
      periapsis = tmp;
    }
    const semiMajorAxis = (apoapsis + periapsis) / 2;
    const eccentricity = apoapsis / semiMajorAxis - 1;
    return new Orbit(referenceBody, semiMajorAxis, eccentricity, inclinationDegs, longitudeOfAscendingNodeDegs, argumentOfPeriapsisDegs, meanAnomalyAtEpoch, timeOfPeriapsisPassage);
  }

  static fromPositionAndVelocity(referenceBody: CelestialBody, position: Vector3, velocity: Vector3, t: number) {
    const mu = referenceBody.gravitationalParameter;
    const r = normV(position);
    const v = normV(velocity);
    const specificAngularMomentum = crossVV(position, velocity);
    let nodeVector: Vector3;
    if (specificAngularMomentum[0] !== 0 || specificAngularMomentum[1] !== 0) {
      nodeVector = normalizeV([-specificAngularMomentum[1], specificAngularMomentum[0], 0]);
    } else {
      nodeVector = [1, 0, 0];
    }
    const eccentricityVector = mulVS(subVV(mulVS(position, v * v - mu / r), mulVS(velocity, dotVV(position, velocity))), 1 / mu);
    const semiMajorAxis = 1 / (2 / r - (v * v) / mu);
    const eccentricity = normV(eccentricityVector);
    const orbit = new Orbit(referenceBody, semiMajorAxis, eccentricity);
    orbit.inclination = Math.acos(specificAngularMomentum[2] / normV(specificAngularMomentum));
    if (eccentricity === 0) {
      orbit.argumentOfPeriapsis = 0;
      orbit.longitudeOfAscendingNode = 0;
    } else {
      orbit.longitudeOfAscendingNode = Math.acos(nodeVector[0]);
      if (nodeVector[1] < 0) {
        orbit.longitudeOfAscendingNode = TWO_PI - orbit.longitudeOfAscendingNode;
      }
      orbit.argumentOfPeriapsis = Math.acos(dotVV(nodeVector, eccentricityVector) / eccentricity);
      if (eccentricityVector[2] < 0) {
        orbit.argumentOfPeriapsis = TWO_PI - orbit.argumentOfPeriapsis;
      }
    }
    let trueAnomaly = Math.acos(dotVV(eccentricityVector, position) / (eccentricity * r));
    if (dotVV(position, velocity) < 0) {
      trueAnomaly = -trueAnomaly;
    }
    const meanAnomaly = orbit.meanAnomalyAtTrueAnomaly(trueAnomaly);
    orbit.timeOfPeriapsisPassage = t - meanAnomaly / orbit.meanMotion();
    return orbit;
  }
}

const circularToEscapeDeltaV = (body: CelestialBody, v0: number, vsoi: number, relativeInclination: number): number => {
  if (!(body instanceof OrbitingCelestialBody)) {
    throw Error("Body has no sphere of influence");
  }
  const mu = body.gravitationalParameter;
  const rsoi = body.sphereOfInfluence;
  const v1 = Math.sqrt(vsoi * vsoi + 2 * v0 * v0 - (2 * mu) / rsoi);
  const r0 = mu / (v0 * v0);
  const e = (r0 * v1 * v1) / mu - 1;
  const ap = (r0 * (1 + e)) / (1 - e);
  if (ap > 0 && ap <= rsoi) {
    return NaN;
  }
  if (relativeInclination) {
    return Math.sqrt(v0 * v0 + v1 * v1 - 2 * v0 * v1 * Math.cos(relativeInclination));
  } else {
    return v1 - v0;
  }
};

const insertionToCircularDeltaV = (body: OrbitingCelestialBody, vsoi: number, v0: number): number => {
  const mu = body.gravitationalParameter;
  const rsoi = body.sphereOfInfluence;
  return Math.sqrt(vsoi * vsoi + 2 * v0 * v0 - (2 * mu) / rsoi) - v0;
};

const ejectionAngleFromPeriapsis = (originBody: OrbitingCelestialBody, initialOrbitalVelocity: number, ejectionDeltaV: number): number => {
  const mu = originBody.gravitationalParameter;
  const rsoi = originBody.sphereOfInfluence;
  const initialOrbitRadius = mu / (initialOrbitalVelocity * initialOrbitalVelocity);
  const v1 = Math.sqrt(ejectionDeltaV * ejectionDeltaV + 2 * initialOrbitalVelocity * initialOrbitalVelocity - (2 * mu) / rsoi);
  const e = (initialOrbitRadius * v1 * v1) / mu - 1;
  const a = initialOrbitRadius / (1 - e);
  const theta = Math.acos((a * (1 - e * e) - rsoi) / (e * rsoi));
  return theta + Math.asin((v1 * initialOrbitRadius) / (ejectionDeltaV * rsoi));
};

const ejectionPeriapsisDirection = (vsoi: Vector3, theta: number): Vector3 => {
  const cosTheta = Math.cos(theta);
  const g = -vsoi[0] / vsoi[1];
  const a = 1 + g * g;
  const b = (2 * g * cosTheta) / vsoi[1];
  const c = (cosTheta * cosTheta) / (vsoi[1] * vsoi[1]) - 1;
  let q: number;
  if (b < 0) {
    q = -0.5 * (b - Math.sqrt(b * b - 4 * a * c));
  } else {
    q = -0.5 * (b + Math.sqrt(b * b - 4 * a * c));
  }
  let vx = q / a;
  let vy = g * vx + cosTheta / vsoi[1];
  if (sign(crossVV([vx, vy, 0], vsoi)[2]) < 0) {
    vx = c / q;
    vy = g * vx + cosTheta / vsoi[1];
  }
  return [vx, vy, 0];
};

const ejectionAngleToPrograde = (periapsis: Vector3, prograde: Vector3): number => {
  prograde = normalizeV([prograde[0], prograde[1], 0]);
  if (crossVV(periapsis, prograde)[2] < 0) {
    return TWO_PI - Math.acos(dotVV(periapsis, prograde));
  } else {
    return Math.acos(dotVV(periapsis, prograde));
  }
};

const ejectionDetails = (ejectionDeltaVector: Vector3, ejectionDeltaV: number, originBody: OrbitingCelestialBody, initialOrbitalVelocity: number, progradeDirection: Vector3) => {
  const ejectionDirection = divVS(ejectionDeltaVector, ejectionDeltaV);
  const theta = ejectionAngleFromPeriapsis(originBody, initialOrbitalVelocity, ejectionDeltaV);
  if (Math.abs(Math.sin(theta)) < Math.abs(ejectionDirection[2])) {
    return [NaN, NaN, NaN];
  } else {
    const periapsisDirection = ejectionPeriapsisDirection(ejectionDirection, theta);
    const ejectionAngle = ejectionAngleToPrograde(periapsisDirection, progradeDirection);
    let ejectionInclination = Math.acos(normalizeV(crossVV(periapsisDirection, ejectionDirection))[2]);
    ejectionInclination *= sign(Math.PI - theta) * sign(ejectionDirection[2]);
    // @todo I think ejectionDeltaV ought to have a different name, because it seems to be being used twice with a different meaning
    ejectionDeltaV = circularToEscapeDeltaV(originBody, initialOrbitalVelocity, ejectionDeltaV, ejectionInclination);
    return [ejectionDeltaV, ejectionInclination, ejectionAngle];
  }
};

enum TransferType {
  BALLISTIC = "ballistic",
  OPTIMAL_PLANE_CHANGE = "optimalPlaneChange",
  OPTIMAL = "optimal",
  PLANE_CHANGE = "planeChange",
}

type Transfer = {
  angle: number;
  orbit?: Orbit;
  ejectionVelocity: Vector3;
  ejectionDeltaVector: Vector3;
  ejectionInclination: number;
  ejectionAngle?: number;
  ejectionDeltaV: number;
  planeChangeAngleToIntercept?: number;
  planeChangeDeltaV: number;
  planeChangeTime?: number;
  planeChangeAngle?: number;
  insertionVelocity: Vector3;
  insertionInclination?: number;
  insertionDeltaV: number;
  deltaV: number;
};

const findTransfer = (transferType: TransferType, opts: TransferOptions): Transfer => {
  switch (transferType) {
    case TransferType.OPTIMAL:
      const ballisticTransfer = findTransfer(TransferType.BALLISTIC, opts);
      if (ballisticTransfer.angle <= HALF_PI) {
        return ballisticTransfer;
      }
      const planeChangeTransfer = findTransfer(TransferType.OPTIMAL_PLANE_CHANGE, opts);
      return ballisticTransfer.deltaV < planeChangeTransfer.deltaV ? ballisticTransfer : planeChangeTransfer;
    case TransferType.OPTIMAL_PLANE_CHANGE:
      const x = getPlaneChangeAngleToIntercept(opts);
      return calculateTransfer(opts, x);
    case TransferType.PLANE_CHANGE:
      return calculateTransfer(opts, HALF_PI);
    case TransferType.BALLISTIC:
      return calculateTransfer(opts, undefined);
  }
};

const calculateTransfer = (opts: TransferOptions, planeChangeAngleToIntercept: number | undefined): Transfer => {
  let planeChangeAngle;
  let planeChangeRotation;
  let p1InOriginPlane;
  // @todo we shouldn't need to set these dummy values - they are bound to be set by the code below
  // it just isn't clear that they always will be
  let ejectionVelocity: Vector3 = [0, 0, 0];
  let insertionVelocity: Vector3 = [0, 0, 0];
  let planeChangeDeltaV;
  let planeChangeTime;
  let ejectionInclination;
  let ejectionAngle;
  let ejectionDeltaV;
  let insertionDeltaV;
  let insertionInclination;
  let orbit;

  if (planeChangeAngleToIntercept != null) {
    const relativeInclination = Math.asin(dotVV(opts.p1, opts.n0) / normV(opts.p1));
    planeChangeAngle = Math.atan2(Math.tan(relativeInclination), Math.sin(planeChangeAngleToIntercept));
  }

  // Assume a counter-clockwise transfer around the +z axis
  let transferAngle = Math.acos(dotVV(opts.p0, opts.p1) / (normV(opts.p0) * normV(opts.p1)));
  if (opts.p0[0] * opts.p1[1] - opts.p0[1] * opts.p1[0] < 0) {
    transferAngle = TWO_PI - transferAngle;
  }

  if (planeChangeAngleToIntercept && planeChangeAngle && transferAngle > HALF_PI) {
    const planeChangeAxis = rotate(quaternionFromAngleAndAxis(-planeChangeAngleToIntercept, opts.n0), projectToPlane(opts.p1, opts.n0));
    planeChangeRotation = quaternionFromAngleAndAxis(planeChangeAngle, planeChangeAxis);
    p1InOriginPlane = rotate(conjugateQ(planeChangeRotation), opts.p1);

    const solutions = solveLambert(opts.referenceBody.gravitationalParameter, opts.p0, p1InOriginPlane, opts.dt);
    ejectionVelocity = solutions[0].ejectionVelocity;
    insertionVelocity = solutions[0].insertionVelocity;
    orbit = Orbit.fromPositionAndVelocity(opts.referenceBody, opts.p0, ejectionVelocity, opts.t0);
    const planeChangeTrueAnomaly = orbit.trueAnomalyAt(opts.t1) - planeChangeAngleToIntercept;
    planeChangeDeltaV = Math.abs(2 * orbit.speedAtTrueAnomaly(planeChangeTrueAnomaly) * Math.sin(planeChangeAngle / 2));
    if (isNaN(planeChangeDeltaV)) {
      planeChangeDeltaV = 0;
    }
    planeChangeTime = orbit.timeAtTrueAnomaly(planeChangeTrueAnomaly, opts.t0);
    insertionVelocity = rotate(planeChangeRotation, insertionVelocity);
  } else {
    const solutions = solveLambert(opts.referenceBody.gravitationalParameter, opts.p0, opts.p1, opts.dt, 10);
    let minDeltaV = Infinity;
    for (let j = 0, len = solutions.length; j < len; j++) {
      const s = solutions[j];
      let dv = normV(subVV(s.ejectionVelocity, opts.v0));
      // @todo original code was this
      // if (typeof finalOrbitVelocity !== "undefined" && finalOrbitVelocity !== null) {
      // but I think this must be a typo, so I changed to finalOrbitalVelocity
      if (opts.finalOrbitalVelocity != null) {
        dv += normV(subVV(s.insertionVelocity, opts.v1));
      }
      if (dv < minDeltaV) {
        minDeltaV = dv;
        ejectionVelocity = s.ejectionVelocity;
        insertionVelocity = s.insertionVelocity;
        transferAngle = s.transferAngle;
      }
    }
    planeChangeDeltaV = 0;
  }

  const ejectionDeltaVector = subVV(ejectionVelocity, opts.v0);
  ejectionDeltaV = normV(ejectionDeltaVector); // This is actually the hyperbolic excess velocity if ejecting from a parking orbit

  if (opts.initialOrbitalVelocity != null) {
    [ejectionDeltaV, ejectionInclination, ejectionAngle] = ejectionDetails(ejectionDeltaVector, ejectionDeltaV, opts.originBody, opts.initialOrbitalVelocity, normalizeV(opts.v0));
  } else {
    ejectionInclination = Math.asin(ejectionDeltaVector[2] / ejectionDeltaV);
  }

  if (opts.finalOrbitalVelocity != null) {
    const insertionDeltaVector = subVV(insertionVelocity, opts.v1);
    insertionDeltaV = normV(insertionDeltaVector); // This is actually the hyperbolic excess velocity if inserting into a parking orbit
    insertionInclination = Math.asin(insertionDeltaVector[2] / insertionDeltaV);
    if (opts.finalOrbitalVelocity) {
      insertionDeltaV = insertionToCircularDeltaV(opts.destinationBody, insertionDeltaV, opts.finalOrbitalVelocity);
    }
  } else {
    insertionDeltaV = 0;
  }

  return {
    angle: transferAngle,
    orbit,
    ejectionVelocity,
    ejectionDeltaVector,
    ejectionInclination,
    ejectionAngle,
    ejectionDeltaV,
    planeChangeAngleToIntercept,
    planeChangeDeltaV,
    planeChangeTime,
    planeChangeAngle: planeChangeTime != null ? planeChangeAngle : 0,
    insertionVelocity,
    insertionInclination,
    insertionDeltaV,
    deltaV: ejectionDeltaV + planeChangeDeltaV + insertionDeltaV,
  };
};

const getPlaneChangeAngleToIntercept = (opts: TransferOptions): number => {
  let x1: number;
  let x2: number;
  if (normV(opts.p0) > normV(opts.p1)) {
    // Transferring to a lower orbit, optimum time to change inclination is 90 degrees to intercept or sooner
    x1 = HALF_PI;
    x2 = Math.PI;
  } else {
    x1 = 0;
    x2 = HALF_PI;
  }
  // This calculates an approximation of the optimal angle to intercept to perform the plane change.
  // The approximation does not take into account the change in the transfer orbit due to the change
  // in the target position rotated into the origin plane as the plane change axis changes.
  // This approximation should be valid so long as the transfer orbit's semi-major axis and eccentricity
  // does not change significantly with the change in the plane change axis.
  const relativeInclination = Math.asin(dotVV(opts.p1, opts.n0) / normV(opts.p1));
  let planeChangeRotation = quaternionFromAngleAndAxis(-relativeInclination, crossVV(opts.p1, opts.n0));
  let p1InOriginPlane = rotate(planeChangeRotation, opts.p1);
  let v1InOriginPlane = rotate(planeChangeRotation, opts.v1);
  let ejectionVelocity = solveLambert(opts.referenceBody.gravitationalParameter, opts.p0, p1InOriginPlane, opts.dt)[0].ejectionVelocity;
  let orbit = Orbit.fromPositionAndVelocity(opts.referenceBody, opts.p0, ejectionVelocity, opts.t0);
  let trueAnomalyAtIntercept = orbit.trueAnomalyAtPosition(p1InOriginPlane);
  let x = goldenSectionSearch(x1, x2, 1e-2, (_x: number) => {
    const _planeChangeAngle = Math.atan2(Math.tan(relativeInclination), Math.sin(_x));
    return Math.abs(2 * orbit.speedAtTrueAnomaly(trueAnomalyAtIntercept - _x) * Math.sin(0.5 * _planeChangeAngle));
  });

  // Refine the initial estimate by running the algorithm again
  const planeChangeAngle = Math.atan2(Math.tan(relativeInclination), Math.sin(x));
  const planeChangeAxis = rotate(quaternionFromAngleAndAxis(-x, opts.n0), projectToPlane(opts.p1, opts.n0));
  planeChangeRotation = quaternionFromAngleAndAxis(planeChangeAngle, planeChangeAxis);
  p1InOriginPlane = rotate(planeChangeRotation, opts.p1);
  v1InOriginPlane = rotate(planeChangeRotation, opts.v1);
  ejectionVelocity = solveLambert(opts.referenceBody.gravitationalParameter, opts.p0, p1InOriginPlane, opts.dt)[0].ejectionVelocity;
  orbit = Orbit.fromPositionAndVelocity(opts.referenceBody, opts.p0, ejectionVelocity, opts.t0);
  trueAnomalyAtIntercept = orbit.trueAnomalyAtPosition(p1InOriginPlane);
  x = goldenSectionSearch(x1, x2, 1e-2, (_x: number) => {
    const _planeChangeAngle = Math.atan2(Math.tan(relativeInclination), Math.sin(_x));
    return Math.abs(2 * orbit.speedAtTrueAnomaly(trueAnomalyAtIntercept - _x) * Math.sin(0.5 * _planeChangeAngle));
  });

  return x;
};

/* @todo come back to this
const transferDetails = (transfer: Transfer, originBody: CelestialBody, t0: number, initialOrbitalVelocity: number) => {
    var burnDirection, ejectionDeltaV, ejectionDeltaVector, ejectionInclination, mu, n0, normalDeltaV, nu0, p0, positionDirection, progradeDeltaV, progradeDirection, radialDeltaV, referenceBody, rsoi, v0, v1, vsoi;
    referenceBody = originBody.orbit.referenceBody;
    nu0 = originBody.orbit.trueAnomalyAt(t0);
    p0 = originBody.orbit.positionAtTrueAnomaly(nu0);
    v0 = originBody.orbit.velocityAtTrueAnomaly(nu0);
    if (transfer.orbit == null) {
        transfer.orbit = Orbit.fromPositionAndVelocity(referenceBody, p0, transfer.ejectionVelocity, t0);
    }
    ejectionDeltaVector = transfer.ejectionDeltaVector;
    ejectionInclination = transfer.ejectionInclination;
    if (initialOrbitalVelocity) {
        mu = originBody.gravitationalParameter;
        rsoi = originBody.sphereOfInfluence;
        vsoi = numeric.normV(ejectionDeltaVector);
        v1 = Math.sqrt(vsoi * vsoi + 2 * initialOrbitalVelocity * initialOrbitalVelocity - 2 * mu / rsoi);
        transfer.ejectionNormalDeltaV = v1 * Math.sin(ejectionInclination);
        transfer.ejectionProgradeDeltaV = v1 * Math.cos(ejectionInclination) - initialOrbitalVelocity;
        transfer.ejectionHeading = Math.atan2(transfer.ejectionProgradeDeltaV, transfer.ejectionNormalDeltaV);
    } else {
        ejectionDeltaV = transfer.ejectionDeltaV;
        positionDirection = normalize(p0);
        progradeDirection = normalize(v0);
        n0 = originBody.orbit.normalVector();
        burnDirection = numeric.divVS(ejectionDeltaVector, ejectionDeltaV);
        transfer.ejectionPitch = Math.asin(numeric.dotVV(burnDirection, positionDirection));
        transfer.ejectionHeading = angleInPlane([0, 0, 1], burnDirection, positionDirection);
        progradeDeltaV = numeric.dotVV(ejectionDeltaVector, progradeDirection);
        normalDeltaV = numeric.dotVV(ejectionDeltaVector, n0);
        radialDeltaV = Math.sqrt(ejectionDeltaV * ejectionDeltaV - progradeDeltaV * progradeDeltaV - normalDeltaV * normalDeltaV);
        if (numeric.dotVV(crossVV(burnDirection, progradeDirection), n0) < 0) {
            radialDeltaV = -radialDeltaV;
        }
        transfer.ejectionProgradeDeltaV = progradeDeltaV;
        transfer.ejectionNormalDeltaV = normalDeltaV;
        transfer.ejectionRadialDeltaV = radialDeltaV;
    }
    return transfer;
};
*/

const refineTransfer = (transfer: Transfer, transferType: TransferType, opts: TransferOptions): Transfer => {
  if (!opts.initialOrbitalVelocity || transfer.ejectionAngle == null) {
    return transfer;
  }
  const originBody = opts.originBody;
  const originOrbit = originBody.orbit;
  const prograde = originOrbit.velocityAtTrueAnomaly(originOrbit.trueAnomalyAt(opts.t0));
  let lastEjectionDeltaVector: Vector3 = transfer.ejectionDeltaVector;

  for (let i = 1; i <= 10; i++) {
    if (isNaN(transfer.deltaV)) {
      return transfer;
    }

    // Calculate the ejection orbit
    const mu = originBody.gravitationalParameter;
    const rsoi = originBody.sphereOfInfluence;
    const vsoi = normV(transfer.ejectionDeltaVector);
    const v0 = Math.sqrt(vsoi * vsoi + 2 * opts.initialOrbitalVelocity * opts.initialOrbitalVelocity - (2 * mu) / rsoi); // Eq 4.15 Velocity at periapsis
    const initialOrbitRadius = mu / (opts.initialOrbitalVelocity * opts.initialOrbitalVelocity);
    const e = (initialOrbitRadius * v0 * v0) / mu - 1; // Ejection orbit eccentricity
    const a = initialOrbitRadius / (1 - e); //  Ejection orbit semi-major axis
    const nu = Math.acos((a * (1 - e * e) - rsoi) / (e * rsoi)); // Eq. 4.82 True anomaly at SOI

    let longitudeOfAscendingNode = Math.atan2(prograde[1], prograde[0]) - transfer.ejectionAngle;
    let argumentOfPeriapsis = 0;
    if (transfer.ejectionInclination < 0) {
      longitudeOfAscendingNode -= Math.PI;
      argumentOfPeriapsis = Math.PI;
    }
    while (longitudeOfAscendingNode < 0) {
      longitudeOfAscendingNode += TWO_PI;
    }

    const ejectionOrbit = new Orbit(opts.originBody, a, e, undefined, undefined, undefined, undefined, opts.t0);
    ejectionOrbit.inclination = transfer.ejectionInclination;
    ejectionOrbit.longitudeOfAscendingNode = longitudeOfAscendingNode;
    ejectionOrbit.argumentOfPeriapsis = argumentOfPeriapsis;

    // Calculate the actual position and time of SoI exit
    const t1 = ejectionOrbit.timeAtTrueAnomaly(nu, opts.t0);
    const originTrueAnomalyAtSOI = originOrbit.trueAnomalyAt(t1);
    const p1 = addVV(ejectionOrbit.positionAtTrueAnomaly(nu), originOrbit.positionAtTrueAnomaly(originTrueAnomalyAtSOI));
    const v1 = originOrbit.velocityAtTrueAnomaly(originTrueAnomalyAtSOI);
    const n1 = normalizeV(crossVV(p1, addVV(v1, transfer.ejectionDeltaVector)));

    // Re-calculate the transfer
    const opts1 = new TransferOptions(opts.originBody, opts.destinationBody, t1, opts.dt - (t1 - opts.t0), 0, opts.finalOrbitalVelocity, p1, v1, n1);
    transfer = findTransfer(transferType, opts1);

    /* tslint:disable-next-line:no-bitwise Bitwise AND is deliberate here as a fast way to test if i is even or odd */
    if (i & 1) {
      lastEjectionDeltaVector = transfer.ejectionDeltaVector;
    } else {
      // Take the average of the last two deltaVectors to avoid diverging series
      transfer.ejectionDeltaVector = mulVS(addVV(lastEjectionDeltaVector, transfer.ejectionDeltaVector), 0.5);
      transfer.ejectionDeltaV = normV(transfer.ejectionDeltaVector);
    }

    // Modify the ejection and total delta-v to take the initial orbit into account
    [transfer.ejectionDeltaV, transfer.ejectionInclination, transfer.ejectionAngle] = ejectionDetails(transfer.ejectionDeltaVector, transfer.ejectionDeltaV, opts.originBody, opts.initialOrbitalVelocity, normalizeV(prograde));
    transfer.deltaV = transfer.ejectionDeltaV + transfer.planeChangeDeltaV + transfer.insertionDeltaV;
  }

  return transfer;
};

const courseCorrection = (transferOrbit: Orbit, destinationOrbit: Orbit, burnTime: number, eta: number) => {
  // Assumes transferOrbit already passes "close" to the destination body at eta
  const mu = transferOrbit.referenceBody.gravitationalParameter;
  const trueAnomaly = transferOrbit.trueAnomalyAt(burnTime);
  const p0 = transferOrbit.positionAtTrueAnomaly(trueAnomaly);
  const v0 = transferOrbit.velocityAtTrueAnomaly(trueAnomaly);
  const n0 = transferOrbit.normalVector();

  const velocityForArrivalAt = (t: number) => {
    const p1 = destinationOrbit.positionAtTrueAnomaly(destinationOrbit.trueAnomalyAt(t));
    return solveLambert(mu, p0, p1, t1 - burnTime)[0].ejectionVelocity;
  };

  // Search for the optimal arrival time within 20% of eta
  const t1Min = Math.max(0.5 * (eta - burnTime), 3600);
  const t1Max = 1.5 * (eta - burnTime);
  let t1 = goldenSectionSearch(t1Min, t1Max, 1e-4, (t: number) => {
    return normSquaredV(subVV(velocityForArrivalAt(burnTime + t), v0));
  });
  t1 = t1 + burnTime; // Convert relative flight time to arrival time

  const correctedVelocity = velocityForArrivalAt(t1);
  const deltaVector = subVV(correctedVelocity, v0);
  const deltaV = normV(deltaVector);
  const burnDirection = divVS(deltaVector, deltaV);
  const positionDirection = normalizeV(p0);
  const pitch = Math.asin(dotVV(burnDirection, positionDirection));
  const heading = angleInPlane([0, 0, 1], burnDirection, positionDirection);
  const progradeDirection = normalizeV(v0);
  const progradeDeltaV = dotVV(deltaVector, progradeDirection);
  const normalDeltaV = dotVV(deltaVector, n0);
  let radialDeltaV = Math.sqrt(deltaV * deltaV - progradeDeltaV * progradeDeltaV - normalDeltaV * normalDeltaV);
  if (dotVV(crossVV(burnDirection, progradeDirection), n0) < 0) {
    radialDeltaV = -radialDeltaV;
  }
  return {
    correctedVelocity,
    deltaVector,
    deltaV,
    pitch,
    heading,
    progradeDeltaV,
    normalDeltaV,
    radialDeltaV,
    arrivalTime: t1,
  };
};

export { Orbit, TransferType, findTransfer };