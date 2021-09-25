import { Orbit } from "./orbit";
import { CelestialBody } from "./celestial-body";
declare class OrbitingCelestialBody extends CelestialBody {
    orbit: Orbit;
    sphereOfInfluence: number;
    constructor(name: string, mass: number, radius: number, siderealRotation: number, orbit: Orbit, atmPressure?: number, atmScaleHeight?: number);
    static fromJSON(json: any): OrbitingCelestialBody;
}
export { OrbitingCelestialBody };
