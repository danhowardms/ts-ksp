import {OrbitingCelestialBodyJSON, makeOrbitingCelestialBody, OrbitingCelestialBody} from "./orbiting-celestial-body";
import {TransferOptions} from "./transfer-options";
import {findTransfer, TransferType, Transfer, Orbit} from "./orbit";
import {
    addVV,
    crossVV,
    dotVV,
    mulVS,
    normalizeV,
    normSquaredV,
    normV,
    signedAngleInPlaneBetween,
    subVV,
    Vector3
} from "./vector3";
import {brentsMethod, goldenSectionSearch} from "./roots";

type SlingshotOptions = {
    startTime: number,
    originBody: OrbitingCelestialBody | OrbitingCelestialBodyJSON,
    slingshotBody: OrbitingCelestialBody | OrbitingCelestialBodyJSON,
    destinationBody: OrbitingCelestialBody | OrbitingCelestialBodyJSON,
    totalDuration: number,
    originOrbitalSpeed: number,
    destinationOrbitalSpeed: number,
};

class SlingshotPlanner {
    startTime: number;
    endTime: number;
    originBody: OrbitingCelestialBody;
    slingshotBody: OrbitingCelestialBody;
    destinationBody: OrbitingCelestialBody;
    originOrbitalSpeed: number;
    destinationOrbitalSpeed: number;
    totalDuration: number;
    firstLegDuration: number = 0;
    secondLegDuration: number = 0;
    slingshotTime: number = 0;

    _transfer1?: Transfer;
    _transfer2?: Transfer;
    _vT1EndRef?: Vector3;
    _vT2StartRef?: Vector3;
    _vSlingshotBodyRef?: Vector3;
    _vT1EndSs?: Vector3;
    _vT2StartSs?: Vector3;
    _ssManeuverPlane?: Vector3;
    _slingshotApproachOrbit?: Orbit;
    _slingshotExitOrbit?: Orbit;
    _pManSs?: Vector3;
    _speedManEndSs?: number;
    _vManStartSs?: Vector3;
    _vManEndSs?: Vector3;
    _deltaVectorMan?: Vector3;
    _deltaVMan?: number;
    _slingshotSoiInTime?: number;
    _slingshotSoiOutTime?: number;

    constructor(opts: SlingshotOptions, ratio?: number) {
        this.startTime = opts.startTime;
        this.endTime = opts.startTime + opts.totalDuration;
        this.originBody = makeOrbitingCelestialBody(opts.originBody);
        this.slingshotBody = makeOrbitingCelestialBody(opts.slingshotBody);
        this.destinationBody = makeOrbitingCelestialBody(opts.destinationBody);
        this.originOrbitalSpeed = opts.destinationOrbitalSpeed;
        this.destinationOrbitalSpeed = opts.destinationOrbitalSpeed;
        this.totalDuration = opts.totalDuration;
        this.setRatio(ratio == null ? 0.5 : ratio);
    }

    reset() {
        this._transfer1 = undefined;
        this._transfer2 = undefined;
        this._vT1EndRef = undefined;
        this._vT2StartRef = undefined;
        this._vSlingshotBodyRef = undefined;
        this._vT1EndSs = undefined;
        this._vT2StartSs = undefined;
        this._ssManeuverPlane = undefined;
        this._slingshotApproachOrbit = undefined;
        this._slingshotExitOrbit = undefined;
        this._pManSs = undefined;
        this._speedManEndSs = undefined;
        this._vManStartSs = undefined;
        this._vManEndSs = undefined;
        this._deltaVectorMan = undefined;
        this._deltaVMan = undefined;
        this._slingshotSoiInTime = undefined;
        this._slingshotSoiOutTime = undefined;
    }

    setRatio(ratio: number) {
        this.reset();
        this.firstLegDuration = this.totalDuration * ratio;
        this.secondLegDuration = this.totalDuration * (1 - ratio);
        this.slingshotTime = this.startTime + this.firstLegDuration;
    }

    // Transfer from originBody to slingshotBody
    get transfer1(): Transfer {
        if (this._transfer1 == null) {
            const transferOpts = new TransferOptions(this.originBody, this.slingshotBody, this.startTime, this.firstLegDuration, this.originOrbitalSpeed, null);
            this._transfer1 = findTransfer(TransferType.OPTIMAL, transferOpts);
        }
        return this._transfer1;
    }

    // Transfer from slingshotBody to destinationBody
    get transfer2(): Transfer {
        if (this._transfer2 == null) {
            const transferOpts = new TransferOptions(this.slingshotBody, this.destinationBody, this.slingshotTime, this.secondLegDuration, null, this.destinationOrbitalSpeed);
            this._transfer2 = findTransfer(TransferType.BALLISTIC, transferOpts);
        }
        return this._transfer2;
    }

    // Velocity vector at end of transfer1, in the system reference frame
    get vT1EndRef(): Vector3 {
        if (this._vT1EndRef == null) {
            const tA1: number = this.transfer1.orbit.trueAnomalyAt(this.slingshotTime);
            this._vT1EndRef = this.transfer1.orbit.velocityAtTrueAnomaly(tA1);
        }
        return this._vT1EndRef;
    }

    // Velocity vector at start of transfer2, in the system reference frame
    get vT2StartRef(): Vector3 {
        if (this._vT2StartRef == null) {
            const tA1: number = this.transfer2.orbit.trueAnomalyAt(this.slingshotTime);
            this._vT2StartRef = this.transfer2.orbit.velocityAtTrueAnomaly(tA1);
        }
        return this._vT2StartRef;
    }

    // Velocity vector of the slingshotBody at time of slingshot
    get vSlingshotBodyRef(): Vector3 {
        if (this._vSlingshotBodyRef == null) {
            const tAs: number = this.slingshotBody.orbit.trueAnomalyAt(this.slingshotTime);
            this._vSlingshotBodyRef = this.slingshotBody.orbit.velocityAtTrueAnomaly(tAs);
        }
        return this._vSlingshotBodyRef;
    }

    // Velocity vector at the end of transfer 1 relative to the slingshotBody
    get vT1EndSs(): Vector3 {
        if (this._vT1EndSs == null) {
            this._vT1EndSs = subVV(this.vT1EndRef, this.vSlingshotBodyRef);
        }
        return this._vT1EndSs;
    }

    // Velocity vector at the start of transfer 2 relative to the slingshotBody
    get vT2StartSs(): Vector3 {
        if (this._vT2StartSs == null) {
            this._vT2StartSs = subVV(this.vT2StartRef, this.vSlingshotBodyRef);
        }
        return this._vT2StartSs;
    }

    // Plane in which out slingshot maneuver will take place
    get ssManeuverPlane(): Vector3 {
        if (this._ssManeuverPlane == null) {
            this._ssManeuverPlane = normalizeV(crossVV(this.vT1EndSs, this.vT2StartSs));
        }
        return this._ssManeuverPlane;
    }

    // Best orbit (around the slingshotBody) for approaching the slingshot maneuver
    get slingshotApproachOrbit(): Orbit {
        if (this._slingshotApproachOrbit != null) {
            return this._slingshotApproachOrbit;
        }
        // The possible approaches are a family of lines parallel to vIn, within the plane
        // We can parameterize these possible approaches by h, the perpendicular distance from the center line
        // at the moment we enter the soi.
        // Direction of h is perpendicular to the plane and vIn
        // Our other basis vector can be a unit vector in the direction opposite to vIn
        const vIn = this.vT1EndSs;
        const vOutTarget = this.vT2StartSs;
        const R = this.slingshotBody.sphereOfInfluence;
        const basis = {
            i: normalizeV(crossVV(vIn, this.ssManeuverPlane)),
            j: normalizeV(mulVS(vIn, -1))
        };
        const approachOrbit = (h: number): Orbit => {
            const pIn = addVV(mulVS(basis.i, h), mulVS(basis.j, Math.sqrt(R * R - h * h)));
            // we're just using this orbit for intermediate calculations
            // so don't need to set the correct time; t=0 will do
            return Orbit.fromPositionAndVelocity(this.slingshotBody, pIn, vIn, 0);
        };

        // Find the best choice of h by comparing the exit velocity for the approach orbit
        // against our desired exit velocity (while also ruling out approaches that are too close to the planet)
        const testFn = (h: number) => {
            const orbit = approachOrbit(h);
            if (orbit.periapsis() < (this.slingshotBody.atmRadius + 2000)) {
                // too close to planet
                return Infinity;
            }
            // Find the true anomaly at the position where this orbit leaves the soi
            const trueAnomalyOut = orbit.trueAnomalyAtRadiusOutbound(R);
            // Now we can calculate the velocity at pOut and see how close it is to our target exit velocity
            const vOut = orbit.velocityAtTrueAnomaly(trueAnomalyOut);
            return normSquaredV(subVV(vOut, vOutTarget));
        };

        // Due to the choice of order of the cross products defining the plane and the i, j vectors
        // possible choices of h will always be 0 < h < R
        // This range puts us the correct side of the slingshot planet so that the planet's gravity
        // pulls us towards the target exit direction.
        // Negative h would put us the wrong side of the planet pulling us away from the target direction.
        const hBest = goldenSectionSearch(0, R, 0.1, testFn);
        const bestOrbit = this._slingshotApproachOrbit = approachOrbit(hBest);
        bestOrbit.timeOfPeriapsisPassage = this.slingshotTime;
        return bestOrbit;
    };

    get slingshotSoiInTime(): number {
        if (this._slingshotSoiInTime != null) {
            return this._slingshotSoiInTime;
        }
        const ta = this.slingshotApproachOrbit.trueAnomalyAtRadiusInbound(this.slingshotBody.sphereOfInfluence);
        return this._slingshotSoiInTime = this.slingshotApproachOrbit.timeAtTrueAnomaly(ta, null);
    }

    // Get position of our maneuver (periapsis of orbitIn), relative to slingshotBody
    get pManSs(): Vector3 {
        if (this._pManSs != null) {
            return this._pManSs;
        }
        return this._pManSs = this.slingshotApproachOrbit.positionAtTrueAnomaly(0);
    }

    // Velocity just prior to maneuver - while on slingshotApproachOrbit, relative to slingshotBody
    get vManStartSs(): Vector3 {
        if (this._vManStartSs != null) {
            return this._vManStartSs;
        }
        return this._vManStartSs = this.slingshotApproachOrbit.velocityAtTrueAnomaly(0);
    }

    // Speed after maneuver, relative to slingshotBody
    get speedManEndSs(): number {
        if (this._speedManEndSs != null) {
            return this._speedManEndSs;
        }
        // Calculate the specific orbital energy of the desired orbit exiting the slingshotBody
        const R = this.slingshotBody.sphereOfInfluence;
        const mu = this.slingshotBody.gravitationalParameter;
        const EOut = 0.5 * normSquaredV(this.vT2StartSs) - (mu / R);

        // From this we can deduce the magnitude of the velocity that we'll need after
        // completing our maneuver
        const rMan = this.slingshotApproachOrbit.periapsis();
        return this._speedManEndSs = Math.sqrt(2 * (EOut + (mu / rMan)));
    }

    // Velocity after maneuver - now on slingshotExitOrbit, relative to slingshotBody
    get vManEndSs(): Vector3 {
        if (this._vManEndSs != null) {
            return this._vManEndSs;
        }
        // Parametrize the possible exit orbits by arg
        // (orbits in our plane which pass through pManSs at speed speedManEndSs)
        const R = this.slingshotBody.sphereOfInfluence;
        const vOutTarget = this.vT2StartSs;
        const basis = {
            i: normalizeV(vOutTarget),
            j: normalizeV(crossVV(this.ssManeuverPlane, vOutTarget))
        };
        const vMan = (arg: number): Vector3 => {
            return addVV(mulVS(basis.i, this.speedManEndSs * Math.cos(arg)), mulVS(basis.j, this.speedManEndSs * Math.sin(arg)));
        };
        const exitOrbit = (arg: number): Orbit => {
            return Orbit.fromPositionAndVelocity(this.slingshotBody, this.pManSs, vMan(arg), 0);
        };

        // For each possible orbit, deduce the exit velocity
        // when the (signed) angle between vOut and vOutTarget reaches zero, they are parallel
        const testFn = (arg: number): number => {
            const orbit = exitOrbit(arg);
            const trueAnomalyOut = orbit.trueAnomalyAtRadiusOutbound(R);
            const vOut = orbit.velocityAtTrueAnomaly(trueAnomalyOut);
            return signedAngleInPlaneBetween(vOut, vOutTarget, this.ssManeuverPlane);
        }

        const relativeAccuracy = 0.001; // @todo what is this?
        const N = 4;
        const a = (Math.PI * 2) / N;
        let bestDeltaVSquared = Infinity;
        let bestVMan: Vector3 = [NaN, NaN, NaN];
        for (let i = 0; i < 4; i++) {
            const arg0 = i * a - Math.PI;
            const argBest = brentsMethod(arg0, arg0 + a, relativeAccuracy, testFn);
            // @todo I don't know why brentsMethod is sometimes giving answers that are not roots
            // but it is, so we need to also check the answer here
            if (isFinite(argBest) && testFn(argBest) < 0.1) {
                const vManMaybe = vMan(argBest);
                const deltaVSquared = normV(subVV(vManMaybe, this.vManStartSs));
                if (deltaVSquared < bestDeltaVSquared) {
                    bestVMan = vManMaybe;
                    bestDeltaVSquared = deltaVSquared;
                }
            }
        }
        return this._vManEndSs = bestVMan;
    };

    // Orbit for exiting the slingshotBody - should be at the starting velocity for transfer2 when we reach the soi
    get slingshotExitOrbit(): Orbit {
        if (this._slingshotExitOrbit != null) {
            return this._slingshotExitOrbit;
        }
        return this._slingshotExitOrbit = Orbit.fromPositionAndVelocity(this.slingshotBody, this.pManSs, this.vManEndSs, this.slingshotTime);
    }

    get slingshotSoiOutTime(): number {
        if (this._slingshotSoiOutTime != null) {
            return this._slingshotSoiOutTime;
        }
        const ta = this.slingshotExitOrbit.trueAnomalyAtRadiusOutbound(this.slingshotBody.sphereOfInfluence);
        return this._slingshotSoiOutTime = this.slingshotExitOrbit.timeAtTrueAnomaly(ta, null);
    }

    // Delta vector for the slingshot maneuver
    get deltaVectorMan(): Vector3 {
        if (this._deltaVectorMan != null) {
            return this._deltaVectorMan;
        }
        return this._deltaVectorMan = subVV(this.vManEndSs, this.vManStartSs);
    }

    get deltaVMan(): number {
        if (this._deltaVMan != null) {
            return this._deltaVMan;
        }
        return this._deltaVMan = normV(this.deltaVectorMan);
    }

    get totalDeltaV(): number {
        return this.transfer1.deltaV + this.deltaVMan + this.transfer2.planeChangeDeltaV + this.transfer2.insertionDeltaV;
    }
}

const findSlingshotRoute = (opts: SlingshotOptions): SlingshotPlanner => {
    const planner = new SlingshotPlanner(opts);
    let bestS = 0.5;
    let bestDeltaV = Infinity;

    const refine = (sMid: number, interval: number, n: number) => {
        const s0 = sMid - 0.5 * interval;
        const ds = interval / n;
        for (let i = 0; i <= n; i++) {
            const s = s0 + i * ds;
            if (s <= 0 || s >= 1) {
                continue;
            }
            planner.setRatio(s);
            let deltaV: number;
            try {
                deltaV = planner.totalDeltaV;
            } catch (e) {
                // Couldn't solve
                deltaV = NaN;
            }
            if (deltaV < bestDeltaV) {
                bestS = s;
                bestDeltaV = planner.totalDeltaV;
            }
        }
    };

    refine(0.5, 0.8, 20);
    refine(bestS, 0.05, 5);

    planner.setRatio(bestS);
    return planner;
};

export {SlingshotOptions, SlingshotPlanner, findSlingshotRoute};