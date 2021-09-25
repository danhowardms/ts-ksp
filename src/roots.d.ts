declare const MACHINE_EPSILON: number;
declare type NumToNumFn = (x: number) => number;
declare const newtonsMethod: (x0: number, f: NumToNumFn, df: NumToNumFn) => number;
declare const brentsMethod: (a: number, b: number, relativeAccuracy: number, f: NumToNumFn, fa?: number | undefined, fb?: number | undefined) => number;
declare const goldenSectionSearch: (x1: number, x2: number, epsilon: number, f: NumToNumFn) => number;
export { MACHINE_EPSILON, newtonsMethod, brentsMethod, goldenSectionSearch };
