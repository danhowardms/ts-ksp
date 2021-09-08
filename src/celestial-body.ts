const G = 6.674e-11;
const TWO_PI = 2 * Math.PI;
const HALF_PI = 0.5 * Math.PI;

const register: Map<string, CelestialBody> = new Map<string, CelestialBody>();

class CelestialBody {
  name: string;
  mass: number;
  radius: number;
  siderealRotation: number;
  atmPressure: number;
  atmScaleHeight: number;
  atmRadius: number;
  gravitationalParameter: number;

  constructor(name: string, mass: number, radius: number, siderealRotation: number, atmPressure?: number, atmScaleHeight?: number) {
    this.name = name;
    this.mass = mass;
    this.radius = radius;
    this.siderealRotation = siderealRotation;
    this.atmPressure = atmPressure == null ? 0 : atmPressure;
    this.atmScaleHeight = atmScaleHeight == null ? 0 : atmScaleHeight;
    this.atmRadius = -Math.log(1e-6) * this.atmScaleHeight + radius;
    this.gravitationalParameter = G * this.mass;

    register.set(name, this);
  }

  circularOrbitVelocity(altitude: number): number {
    return Math.sqrt(this.gravitationalParameter / (altitude + this.radius));
  }

  siderealTimeAt(longitude: number, time: number): number {
    const result = ((time / this.siderealRotation) * TWO_PI + HALF_PI + longitude) % TWO_PI;
    if (result < 0) {
      return result + TWO_PI;
    } else {
      return result;
    }
  }

  /*
  // @todo the children are actually all OrbitingCelestialBody instances, but I can't import OrbitingCelestialBody
  // in this file without creating a circular reference - must be some way round this...
  children(): OrbitingCelestialBody[] {
    const _children: OrbitingCelestialBody[] = [];
    for (const body of register.values()) {
      if (body instanceof OrbitingCelestialBody && body.orbit.referenceBody === this) {
        _children.push(body);
      }
    }
    return _children;
  }
  */

  static getByName(name: string): CelestialBody | undefined {
    return register.get(name);
  }
}

export { CelestialBody };