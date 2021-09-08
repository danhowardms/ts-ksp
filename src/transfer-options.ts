import { CelestialBody } from "./celestial-body";
import { Orbit } from "./orbit";
import { Vector3 } from "./vector3";
import { OrbitingCelestialBody } from "./orbiting-celestial-body";

class TransferOptions {
  referenceBody: CelestialBody;
  originBody: OrbitingCelestialBody;
  destinationBody: OrbitingCelestialBody;
  t0: number;
  dt: number;
  t1: number;
  initialOrbitalVelocity: number | null;
  finalOrbitalVelocity: number | null;
  p0: Vector3;
  v0: Vector3;
  n0: Vector3;
  p1: Vector3;
  v1: Vector3;

  constructor(originBody: OrbitingCelestialBody, destinationBody: OrbitingCelestialBody, t0: number, dt: number, initialOrbitalVelocity: number | null, finalOrbitalVelocity: number | null, p0?: Vector3, v0?: Vector3, n0?: Vector3, p1?: Vector3, v1?: Vector3) {
    this.referenceBody = originBody.orbit.referenceBody;
    this.originBody = originBody;
    this.destinationBody = destinationBody;
    this.t0 = t0;
    this.dt = dt;
    this.t1 = t0 + dt;
    this.initialOrbitalVelocity = initialOrbitalVelocity;
    this.finalOrbitalVelocity = finalOrbitalVelocity;
    [this.p0, this.v0] = TransferOptions.fillInMissing(originBody.orbit, t0, p0, v0);
    [this.p1, this.v1] = TransferOptions.fillInMissing(originBody.orbit, t0 + dt, p1, v1);
    this.n0 = n0 == null ? originBody.orbit.normalVector() : n0;
  }

  static fillInMissing(orbit: Orbit, t: number, p: Vector3 | undefined, v: Vector3 | undefined): [Vector3, Vector3] {
    if (p && v) {
      return [p, v];
    }
    const nu = orbit.trueAnomalyAt(t);
    return [p == null ? orbit.positionAtTrueAnomaly(nu) : p, v == null ? orbit.velocityAtTrueAnomaly(nu) : v];
  }
}

export { TransferOptions };
