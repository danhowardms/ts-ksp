declare class CelestialBody {
    name: string;
    mass: number;
    radius: number;
    siderealRotation: number;
    atmPressure: number;
    atmScaleHeight: number;
    atmRadius: number;
    gravitationalParameter: number;
    constructor(name: string, mass: number, radius: number, siderealRotation: number, atmPressure?: number, atmScaleHeight?: number);
    /**
     * Calculate the speed required to maintain a circular orbit at the given altitude above the body's surface
     *
     * @param altitude
     */
    circularOrbitVelocity(altitude: number): number;
    siderealTimeAt(longitude: number, time: number): number;
    static getByName(name: string): CelestialBody | undefined;
    static fromJSON(json: any): CelestialBody;
}
export { CelestialBody };
