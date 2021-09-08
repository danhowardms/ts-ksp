type Opaque<T, K> = T & { __opaque__: K };
type AngleDegrees = Opaque<number, "AngleDegrees">;

export { AngleDegrees };
