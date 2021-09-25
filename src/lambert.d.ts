import { Vector3 } from "./vector3";
declare type LambertSolution = {
    ejectionVelocity: Vector3;
    insertionVelocity: Vector3;
    transferAngle: number;
};
/**
 * Find solutions to Lambert's problem:
 * Find the orbit that passes through two given points at two given times
 *
 * @param mu Rate of acceleration towards the central body (m/s2 I think)
 * @param pos1 Position vector at the start of the transfer
 * @param pos2 Position vector at the end of the transfer
 * @param dt Time to move from pos1 to pos2
 * @param maxRevs Max number of revolutions around an elliptical transfer orbit to consider
 * @param prograde ??? I'm not sure what this parameter means exactly, but it's always omitted where it was used in
 *        AlexMoon's transfer calculator, which means it would always take the default value, 1. Maybe -1 would give
 *        us a retrograde transfer?
 */
declare const solveLambert: (mu: number, pos1: Vector3, pos2: Vector3, dt: number, maxRevs?: number | undefined, prograde?: number | undefined) => LambertSolution[];
export { LambertSolution, solveLambert };
