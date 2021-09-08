import { CelestialBody } from "./celestial-body";
import { Orbit } from "./orbit";
import { AngleDegrees } from "./utility-types";
import { OrbitingCelestialBody } from "./orbiting-celestial-body";

const kerbol = new CelestialBody("Kerbol", 1.756567e28, 2.616e8, 432000);
const moho = new OrbitingCelestialBody("Moho", 2.5263617e21, 250000, 1210000, new Orbit(kerbol, 5263138304, 0.2, 7.0 as AngleDegrees, 70.0 as AngleDegrees, 15.0 as AngleDegrees, 3.14));
const eve = new OrbitingCelestialBody("Eve", 1.2244127e23, 700000, 80500, new Orbit(kerbol, 9832684544, 0.01, 2.1 as AngleDegrees, 15.0 as AngleDegrees, 0 as AngleDegrees, 3.14), 5, 7000);
const gilly = new OrbitingCelestialBody("Gilly", 1.2420512e17, 13000, 28255, new Orbit(eve, 31500000, 0.55, 12.0 as AngleDegrees, 80.0 as AngleDegrees, 10.0 as AngleDegrees, 0.9));
const kerbin = new OrbitingCelestialBody("Kerbin", 5.2915793e22, 600000, 21600, new Orbit(kerbol, 13599840256, 0.0, 0 as AngleDegrees, 0 as AngleDegrees, 0 as AngleDegrees, 3.14), 1, 5000);
const mun = new OrbitingCelestialBody("Mun", 9.7600236e20, 200000, 138984.38, new Orbit(kerbin, 12000000, 0.0, 0 as AngleDegrees, 0 as AngleDegrees, 0 as AngleDegrees, 1.7));
const minmus = new OrbitingCelestialBody("Minmus", 2.6457897e19, 60000, 40400, new Orbit(kerbin, 47000000, 0.0, 6.0 as AngleDegrees, 78.0 as AngleDegrees, 38.0 as AngleDegrees, 0.9));
const duna = new OrbitingCelestialBody("Duna", 4.5154812e21, 320000, 65517.859, new Orbit(kerbol, 20726155264, 0.051, 0.06 as AngleDegrees, 135.5 as AngleDegrees, 0 as AngleDegrees, 3.14), 0.2, 3000);
const ike = new OrbitingCelestialBody("Ike", 2.7821949e20, 130000, 65517.862, new Orbit(duna, 3200000, 0.03, 0.2 as AngleDegrees, 0 as AngleDegrees, 0 as AngleDegrees, 1.7));
const dres = new OrbitingCelestialBody("Dres", 3.2191322e20, 138000, 34800, new Orbit(kerbol, 40839348203, 0.145, 5.0 as AngleDegrees, 280.0 as AngleDegrees, 90.0 as AngleDegrees, 3.14));
const jool = new OrbitingCelestialBody("Jool", 4.2332635e24, 6000000, 36000, new Orbit(kerbol, 68773560320, 0.05, 1.304 as AngleDegrees, 52.0 as AngleDegrees, 0 as AngleDegrees, 0.1), 15, 10000);
const laythe = new OrbitingCelestialBody("Laythe", 2.9397663e22, 500000, 52980.879, new Orbit(jool, 27184000, 0, 0 as AngleDegrees, 0 as AngleDegrees, 0 as AngleDegrees, 3.14), 0.8, 4000);
const vall = new OrbitingCelestialBody("Vall", 3.1088028e21, 300000, 105962.09, new Orbit(jool, 43152000, 0, 0 as AngleDegrees, 0 as AngleDegrees, 0 as AngleDegrees, 0.9));
const tylo = new OrbitingCelestialBody("Tylo", 4.2332635e22, 600000, 211926.36, new Orbit(jool, 68500000, 0, 0.025 as AngleDegrees, 0 as AngleDegrees, 0 as AngleDegrees, 3.14));
const bop = new OrbitingCelestialBody("Bop", 3.7261536e19, 65000, 544507.4, new Orbit(jool, 128500000, 0.235, 15.0 as AngleDegrees, 10.0 as AngleDegrees, 25.0 as AngleDegrees, 0.9));
const pol = new OrbitingCelestialBody("Pol", 1.0813636e19, 44000, 901902.62, new Orbit(jool, 179890000, 0.17085, 4.25 as AngleDegrees, 2.0 as AngleDegrees, 15.0 as AngleDegrees, 0.9));
const eeloo = new OrbitingCelestialBody("Eeloo", 1.1149358e21, 210000, 19460, new Orbit(kerbol, 90118820000, 0.26, 6.15 as AngleDegrees, 50.0 as AngleDegrees, 260.0 as AngleDegrees, 3.14));

export { kerbol, moho, eve, gilly, kerbin, mun, minmus, duna, ike, dres, jool, laythe, vall, tylo, bop, pol, eeloo };
