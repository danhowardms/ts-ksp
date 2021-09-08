import {kerbin, moho} from "./kerbol-system";
import {TransferOptions} from "./transfer-options";
import {OrbitingCelestialBody} from "./orbiting-celestial-body";
import {Vector3} from "./vector3";
import {TransferType, findTransfer} from "./orbit";

const opts = new TransferOptions(kerbin, moho, 0, 426 * 6 * 60 * 60, 2600, 800);
const transfer = findTransfer(TransferType.OPTIMAL_PLANE_CHANGE, opts);
console.log(transfer);

