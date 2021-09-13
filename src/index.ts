import {CelestialBody} from "./celestial-body";
import * as KerbolSystem from "./kerbol-system";
import {LambertSolution, solveLambert} from "./lambert";
import {Orbit, TransferType, findTransfer} from "./orbit";
import {OrbitingCelestialBody} from "./orbiting-celestial-body";
import { Quaternion, addQQ, conjugateQ, normalizeQ, concatQQ, quaternionFromAngleAndAxis, quaternionFromStartAndEndVectors, vectorToQuaternion, quaternionToVector, rotate } from './quaternion';
import {newtonsMethod, brentsMethod, goldenSectionSearch} from "./roots";
import {SlingshotOptions, findSlingshotRoute} from "./slingshot-planner";
import {TransferOptions} from "./transfer-options";
import {AngleDegrees} from "./utility-types";
import {Vector3, addVV, subVV, mulVS, divVS, normSquaredV, normV, normalizeV, dotVV, crossVV, projectToPlane, angleBetween, signedAngleInPlaneBetween} from "./vector3";

export {
    CelestialBody,
    OrbitingCelestialBody,
    KerbolSystem,
    LambertSolution,
    Orbit, TransferType, findTransfer,
    Quaternion, addQQ, conjugateQ, normalizeQ, concatQQ, quaternionFromAngleAndAxis, quaternionFromStartAndEndVectors, vectorToQuaternion, quaternionToVector, rotate,
    newtonsMethod, brentsMethod, goldenSectionSearch,
    SlingshotOptions, findSlingshotRoute,
    TransferOptions,
    AngleDegrees,
    Vector3, addVV, subVV, mulVS, divVS, normSquaredV, normV, normalizeV, dotVV, crossVV, projectToPlane, angleBetween, signedAngleInPlaneBetween
};
