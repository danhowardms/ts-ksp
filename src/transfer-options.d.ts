import { CelestialBody } from "./celestial-body";
import { Orbit } from "./orbit";
import { Vector3 } from "./vector3";
import { OrbitingCelestialBody } from "./orbiting-celestial-body";
declare class TransferOptions {
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
    constructor(originBody: OrbitingCelestialBody, destinationBody: OrbitingCelestialBody, t0: number, dt: number, initialOrbitalVelocity: number | null, finalOrbitalVelocity: number | null, p0?: Vector3, v0?: Vector3, n0?: Vector3, p1?: Vector3, v1?: Vector3);
    static fillInMissing(orbit: Orbit, t: number, p: Vector3 | undefined, v: Vector3 | undefined): [Vector3, Vector3];
}
export { TransferOptions };
