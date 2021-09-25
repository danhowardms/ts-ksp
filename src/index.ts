import {CelestialBodyJSON, serializeCelestialBody, makeCelestialBody, CelestialBody} from "./celestial-body";
import * as KerbolSystem from "./kerbol-system";
import {LambertSolution, solveLambert} from "./lambert";
import {OrbitJSON, serializeOrbit, makeOrbit, Orbit, TransferType, findTransfer} from "./orbit";
import {OrbitingCelestialBodyJSON, serializeOrbitingCelestialBody, makeOrbitingCelestialBody, OrbitingCelestialBody} from "./orbiting-celestial-body";
import { Quaternion, addQQ, conjugateQ, normalizeQ, concatQQ, quaternionFromAngleAndAxis, quaternionFromStartAndEndVectors, vectorToQuaternion, quaternionToVector, rotate } from './quaternion';
import {newtonsMethod, brentsMethod, goldenSectionSearch} from "./roots";
import {SlingshotOptions, findSlingshotRoute, SlingshotPlanner} from "./slingshot-planner";
import {TransferOptions} from "./transfer-options";
import {AngleDegrees} from "./utility-types";
import {Vector3, addVV, subVV, mulVS, divVS, normSquaredV, normV, normalizeV, dotVV, crossVV, projectToPlane, angleBetween, signedAngleInPlaneBetween} from "./vector3";

export {
    CelestialBodyJSON, serializeCelestialBody, makeCelestialBody, CelestialBody,
    OrbitingCelestialBodyJSON, serializeOrbitingCelestialBody, makeOrbitingCelestialBody, OrbitingCelestialBody,
    KerbolSystem,
    LambertSolution,
    OrbitJSON, serializeOrbit, makeOrbit, Orbit, TransferType, findTransfer,
    Quaternion, addQQ, conjugateQ, normalizeQ, concatQQ, quaternionFromAngleAndAxis, quaternionFromStartAndEndVectors, vectorToQuaternion, quaternionToVector, rotate,
    newtonsMethod, brentsMethod, goldenSectionSearch,
    SlingshotOptions, findSlingshotRoute, SlingshotPlanner,
    TransferOptions,
    AngleDegrees,
    Vector3, addVV, subVV, mulVS, divVS, normSquaredV, normV, normalizeV, dotVV, crossVV, projectToPlane, angleBetween, signedAngleInPlaneBetween
};
