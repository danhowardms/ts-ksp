import { OrbitingCelestialBody } from "./orbiting-celestial-body";
import { Transfer, Orbit } from "./orbit";
import { Vector3 } from "./vector3";
declare type SlingshotOptions = {
    startTime: number;
    originBody: OrbitingCelestialBody;
    slingshotBody: OrbitingCelestialBody;
    destinationBody: OrbitingCelestialBody;
    totalDuration: number;
    originOrbitalSpeed: number;
    destinationOrbitalSpeed: number;
};
declare class SlingshotPlanner {
    startTime: number;
    endTime: number;
    originBody: OrbitingCelestialBody;
    slingshotBody: OrbitingCelestialBody;
    destinationBody: OrbitingCelestialBody;
    originOrbitalSpeed: number;
    destinationOrbitalSpeed: number;
    totalDuration: number;
    firstLegDuration: number;
    secondLegDuration: number;
    slingshotTime: number;
    _transfer1?: Transfer;
    _transfer2?: Transfer;
    _vT1EndRef?: Vector3;
    _vT2StartRef?: Vector3;
    _vSlingshotBodyRef?: Vector3;
    _vT1EndSs?: Vector3;
    _vT2StartSs?: Vector3;
    _ssManeuverPlane?: Vector3;
    _slingshotApproachOrbit?: Orbit;
    _slingshotExitOrbit?: Orbit;
    _pManSs?: Vector3;
    _speedManEndSs?: number;
    _vManStartSs?: Vector3;
    _vManEndSs?: Vector3;
    _deltaVectorMan?: Vector3;
    _deltaVMan?: number;
    _slingshotSoiInTime?: number;
    _slingshotSoiOutTime?: number;
    constructor(opts: SlingshotOptions, ratio?: number);
    reset(): void;
    setRatio(ratio: number): void;
    get transfer1(): Transfer;
    get transfer2(): Transfer;
    get vT1EndRef(): Vector3;
    get vT2StartRef(): Vector3;
    get vSlingshotBodyRef(): Vector3;
    get vT1EndSs(): Vector3;
    get vT2StartSs(): Vector3;
    get ssManeuverPlane(): Vector3;
    get slingshotApproachOrbit(): Orbit;
    get slingshotSoiInTime(): number;
    get pManSs(): Vector3;
    get vManStartSs(): Vector3;
    get speedManEndSs(): number;
    get vManEndSs(): Vector3;
    get slingshotExitOrbit(): Orbit;
    get slingshotSoiOutTime(): number;
    get deltaVectorMan(): Vector3;
    get deltaVMan(): number;
    get totalDeltaV(): number;
}
declare const findSlingshotRoute: (opts: SlingshotOptions) => SlingshotPlanner;
export { SlingshotOptions, SlingshotPlanner, findSlingshotRoute };