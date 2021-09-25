import { Vector3 } from "./vector3";
import { Quaternion } from "./quaternion";
import { TransferOptions } from "./transfer-options";
import { CelestialBody } from "./celestial-body";
import { AngleDegrees } from "./utility-types";
declare class Orbit {
    referenceBody: CelestialBody;
    semiMajorAxis: number;
    eccentricity: number;
    inclination: number;
    longitudeOfAscendingNode: number;
    argumentOfPeriapsis: number;
    meanAnomalyAtEpoch: number;
    timeOfPeriapsisPassage: number | undefined;
    constructor(referenceBody1: CelestialBody, semiMajorAxis1: number, eccentricity1: number, inclinationDegs?: AngleDegrees, longitudeOfAscendingNodeDegs?: AngleDegrees, argumentOfPeriapsisDegs?: AngleDegrees, meanAnomalyAtEpoch1?: number, timeOfPeriapsisPassage1?: number);
    inclinationDegs(): AngleDegrees;
    longitudeOfAscendingNodeDegs(): AngleDegrees;
    argumentOfPeriapsisDegs(): AngleDegrees;
    isHyperbolic(): boolean;
    apoapsis(): number;
    periapsis(): number;
    apoapsisAltitude(): number;
    periapsisAltitude(): number;
    semiMinorAxis(): number;
    semiLatusRectum(): number;
    meanMotion(): number;
    period(): number;
    rotationToReferenceFrame(): Quaternion;
    normalVector(): Vector3;
    phaseAngle(orbit: Orbit, t: number): number;
    meanAnomalyAt(t: number): number;
    eccentricAnomalyAt(t: number): number;
    trueAnomalyAt(t: number): number;
    trueAnomalyAtRadiusOutbound(r: number): number;
    trueAnomalyAtRadiusInbound(r: number): number;
    positionAt(t: number): Vector3;
    eccentricAnomalyAtTrueAnomaly(tA: number): number;
    meanAnomalyAtTrueAnomaly(tA: number): number;
    timeAtTrueAnomaly(tA: number, t0: number | null): number;
    radiusAtTrueAnomaly(tA: number): number;
    altitudeAtTrueAnomaly(tA: number): number;
    speedAtTrueAnomaly(tA: number): number;
    positionAtTrueAnomaly(tA: number): Vector3;
    velocityAtTrueAnomaly(tA: number): Vector3;
    trueAnomalyAtPosition(p: Vector3): number;
    static fromJSON(json: any): Orbit;
    static fromApoapsisAndPeriapsis(referenceBody: CelestialBody, apoapsis: number, periapsis: number, inclinationDegs: AngleDegrees, longitudeOfAscendingNodeDegs: AngleDegrees, argumentOfPeriapsisDegs: AngleDegrees, meanAnomalyAtEpoch: number, timeOfPeriapsisPassage: number): Orbit;
    static fromPositionAndVelocity(referenceBody: CelestialBody, position: Vector3, velocity: Vector3, t: number): Orbit;
}
declare enum TransferType {
    BALLISTIC = "ballistic",
    OPTIMAL_PLANE_CHANGE = "optimalPlaneChange",
    OPTIMAL = "optimal",
    PLANE_CHANGE = "planeChange"
}
declare type Transfer = {
    t0: number;
    t1: number;
    angle: number;
    orbit: Orbit;
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
declare const findTransfer: (transferType: TransferType, opts: TransferOptions) => Transfer;
export { Orbit, TransferType, Transfer, findTransfer };
