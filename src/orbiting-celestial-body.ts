import { Orbit } from "./orbit";
import { CelestialBody } from "./celestial-body";

class OrbitingCelestialBody extends CelestialBody {
  orbit: Orbit;
  sphereOfInfluence: number;

  constructor(name: string, mass: number, radius: number, siderealRotation: number, orbit: Orbit, atmPressure?: number, atmScaleHeight?: number) {
    super(name, mass, radius, siderealRotation, atmPressure, atmScaleHeight);
    this.orbit = orbit;
    this.sphereOfInfluence = orbit.semiMajorAxis * Math.pow(mass / orbit.referenceBody.mass, 0.4);
  }
}

export { OrbitingCelestialBody };
