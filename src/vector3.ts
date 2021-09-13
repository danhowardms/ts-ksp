type Vector3 = [number, number, number];

const addVV = (a: Vector3, b: Vector3): Vector3 => {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
};

const subVV = (a: Vector3, b: Vector3): Vector3 => {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
};

const mulVS = (v: Vector3, s: number): Vector3 => {
  return [s * v[0], s * v[1], s * v[2]];
};

const divVS = (v: Vector3, s: number): Vector3 => {
  return mulVS(v, 1 / s);
};

const normSquaredV = (v: Vector3): number => {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
};

const normV = (v: Vector3) => {
  return Math.sqrt(normSquaredV(v));
};

const normalizeV = (v: Vector3): Vector3 => {
  return divVS(v, normV(v));
};

const dotVV = (a: Vector3, b: Vector3): number => {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
};

const crossVV = (a: Vector3, b: Vector3): Vector3 => {
  const r: Vector3 = [0, 0, 0];
  r[0] = a[1] * b[2] - a[2] * b[1];
  r[1] = a[2] * b[0] - a[0] * b[2];
  r[2] = a[0] * b[1] - a[1] * b[0];
  return r;
};

const angleBetween = (p0: Vector3, p1: Vector3) => {
  return Math.acos(dotVV(p0, p1) / (normV(p0) * normV(p1)));
};

const signedAngleInPlaneBetween = (p0: Vector3, p1: Vector3, planeNormal: Vector3) => {
  return Math.atan2(dotVV(crossVV(p0, p1), planeNormal), dotVV(p0, p1));
};

// Find the projection of the point p onto the plane (through the origin) defined by the unit normal vector n
const projectToPlane = (p: Vector3, n: Vector3): Vector3 => {
  return subVV(p, mulVS(n, dotVV(p, n)));
};

export { Vector3, addVV, subVV, mulVS, divVS, normSquaredV, normV, normalizeV, dotVV, crossVV, angleBetween, signedAngleInPlaneBetween, projectToPlane };
