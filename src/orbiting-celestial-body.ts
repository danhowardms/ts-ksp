import { OrbitJSON, Orbit, serializeOrbit, makeOrbit } from "./orbit";
import { CelestialBody } from "./celestial-body";

type OrbitingCelestialBodyJSON = {
  name: string,
  mass: number,
  radius: number,
  siderealRotation: number,
  orbit: OrbitJSON,
  atmPressure?: number,
  atmScaleHeight?: number,
};

class OrbitingCelestialBody extends CelestialBody {
  orbit: Orbit;
  sphereOfInfluence: number;

  constructor(name: string, mass: number, radius: number, siderealRotation: number, orbit: Orbit, atmPressure?: number, atmScaleHeight?: number) {
    super(name, mass, radius, siderealRotation, atmPressure, atmScaleHeight);
    this.orbit = orbit;
    this.sphereOfInfluence = orbit.semiMajorAxis * Math.pow(mass / orbit.referenceBody.mass, 0.4);
  }
}

const serializeOrbitingCelestialBody = (input: OrbitingCelestialBody | OrbitingCelestialBodyJSON): OrbitingCelestialBodyJSON => {
  if (input instanceof OrbitingCelestialBody) {
    return {
      name: input.name,
      mass: input.mass,
      radius: input.radius,
      siderealRotation: input.siderealRotation,
      orbit: serializeOrbit(input.orbit),
      atmPressure: input.atmPressure,
      atmScaleHeight: input.atmScaleHeight,
    };
  } else {
    return input;
  }
};

const makeOrbitingCelestialBody = (input: OrbitingCelestialBody | OrbitingCelestialBodyJSON): OrbitingCelestialBody => {
  if (input instanceof  OrbitingCelestialBody) {
    return input;
  } else {
    const orbit = makeOrbit(input.orbit);
    return new OrbitingCelestialBody(input.name, input.mass, input.radius, input.siderealRotation, orbit, input.atmPressure, input.atmScaleHeight);
  }
};

export { OrbitingCelestialBodyJSON, serializeOrbitingCelestialBody, makeOrbitingCelestialBody, OrbitingCelestialBody };
