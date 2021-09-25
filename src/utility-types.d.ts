declare type Opaque<T, K> = T & {
    __opaque__: K;
};
declare type AngleDegrees = Opaque<number, "AngleDegrees">;
export { AngleDegrees };
