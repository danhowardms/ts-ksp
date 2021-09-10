import { Vector3, normalizeV, dotVV } from "./vector3";

type Quaternion = [number, number, number, number];

const addQQ = (a: Quaternion, b: Quaternion): Quaternion => {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]];
};

const conjugateQ = (q: Quaternion): Quaternion => {
  return [-q[0], -q[1], -q[2], q[3]];
};

const normalizeQ = (q: Quaternion): Quaternion => {
  const s = 1 / (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  return [s * q[0], s * q[1], s * q[2], s * q[3]];
};

const concatQQ = (q0: Quaternion, q1: Quaternion): Quaternion => {
  const x = q0[3] * q1[0] + q0[0] * q1[3] + q0[1] * q1[2] - q0[2] * q1[1];
  const y = q0[3] * q1[1] - q0[0] * q1[2] + q0[1] * q1[3] + q0[2] * q1[0];
  const z = q0[3] * q1[2] + q0[0] * q1[1] - q0[1] * q1[0] + q0[2] * q1[3];
  const w = q0[3] * q1[3] - q0[0] * q1[0] - q0[1] * q1[1] - q0[2] * q1[2];
  return [x, y, z, w];
};

const quaternionFromAngleAndAxis = (angle: number, axis: Vector3): Quaternion => {
  const halfAngle = 0.5 * angle;
  const sin = Math.sin(halfAngle);
  axis = normalizeV(axis);
  return normalizeQ([sin * axis[0], sin * axis[1], sin * axis[2], Math.cos(halfAngle)]);
};

const quaternionFromStartAndEndVectors = (start: Vector3, end: Vector3): Quaternion => {
  start = normalizeV(start);
  end = normalizeV(end);
  const dot = dotVV(start, end);
  if (dot > 1.0 - 1e-6) {
    return [0, 0, 0, 1];
  } else if (dot < -(1.0 - 1e-6)) {
    return quaternionFromAngleAndAxis(Math.PI, [0, 0, 1]);
  } else {
    const s = Math.sqrt(2 * (1 + dot));
    const invs = 1 / s;
    return normalizeQ([(start[1] * end[2] - start[2] * end[1]) * invs, (start[2] * end[0] - start[0] * end[2]) * invs, (start[0] * end[1] - start[1] * end[0]) * invs, 0.5 * s]);
  }
};

const vectorToQuaternion = (v: Vector3): Quaternion => {
  return [v[0], v[1], v[2], 0];
};

const quaternionToVector = (q: Quaternion): Vector3 => {
  return [q[0], q[1], q[2]];
};

const rotate = (q: Quaternion, v: Vector3): Vector3 => {
  const p = vectorToQuaternion(v);
  return quaternionToVector(concatQQ(concatQQ(q, p), conjugateQ(q)));
};

const quaternionToRotationMatrix = (q: Quaternion) => {
  // First row of the rotation matrix
  const r00 = 2 * (q[0] * q[0] + q[1] * q[1]) - 1;
  const r01 = 2 * (q[1] * q[2] - q[0] * q[3]);
  const r02 = 2 * (q[1] * q[3] + q[0] * q[2]);
  // Second row of the rotation matrix
  const r10 = 2 * (q[1] * q[2] + q[0] * q[3]);
  const r11 = 2 * (q[0] * q[0] + q[2] * q[2]) - 1;
  const r12 = 2 * (q[2] * q[3] - q[0] * q[1]);
  // Third row of the rotation matrix
  const r20 = 2 * (q[1] * q[3] - q[0] * q[2]);
  const r21 = 2 * (q[2] * q[3] + q[0] * q[1]);
  const r22 = 2 * (q[0] * q[0] + q[3] * q[3]) - 1;
  // 3x3 rotation matrix
  return [[r00, r01, r02], [r10, r11, r12], [r20, r21, r22]];
};

export { Quaternion, addQQ, conjugateQ, normalizeQ, concatQQ, quaternionFromAngleAndAxis, quaternionFromStartAndEndVectors, vectorToQuaternion, quaternionToVector, rotate };
