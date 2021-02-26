extern crate ndarray;
extern crate num;
extern crate serde;
extern crate serde_yaml;
extern crate image;

use std::env;
use std::f64::consts::PI;
use std::fs;
use num::traits::{
    Float,
    FromPrimitive,
};

use ndarray::prelude::*;
use ndarray::{
    Array2,
};

use serde::{
    Serialize,
    Deserialize,
};

const DTOR: f64 = PI / 180.0;

// inverse gamma correction
const R0: f64 = 0.0031308;
const B0: f64 = 0.055;
const PWR: f64 = 0.416666;

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct Config {
    input: Input,
    waves: Vec<Wave>,
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct Input {
  qsol: f64,
  psfc: f64,
  galb: f64,
  atau0: f64,
  angexp: f64,
  aomg: f64,
  aasy: f64,
  the0: f64,
  phi0: f64,
  nu: usize,
  nv: usize,
  cphi: f64,
  cthe: f64,
  cpsi: f64,
  caov: f64,
  cpar: f64,
  chgt: f64,
  cimax: f64,
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct Wave {
    min: f64,
    max: f64,
    fsol: f64,
    cmf: Vec3,
}

fn parse_yaml(path: &str) -> Result<Config, serde_yaml::Error> {
    let config_file = fs::read_to_string(path).unwrap();
    serde_yaml::from_str(&config_file)
}

fn validate(config: &Config) -> Result<(), &str> {
    if config.input.the0 < 0.0 || config.input.the0 >= 90.0 {
        return Err("the0 must be in [0,90).")
    }

    Ok(())
}

fn new_with_cartesian(r: f64, the: f64, phi: f64) -> Array1<f64> {
    let cosq = the.cos();
    let sinq = the.sin();
    arr1(&[
         r * sinq * phi.cos(),
         r * sinq * phi.sin(),
         r * cosq,
    ])
}

fn rotate_matrix_y(ang: f64) -> Array2<f64> {
    let z = 0.0;
    let u = 1.0;
    let s = ang.sin();
    let c = ang.cos();
    arr2(&[
         [c, z, s],
         [z, u, z],
         [-s, z, c],
    ])
}

fn rotate_matrix_z(ang: f64) -> Array2<f64> {
    let z = 0.0;
    let u = 1.0;
    let s = ang.sin();
    let c = ang.cos();
    arr2(&[
         [c, -s, z],
         [s, c, z],
         [z, z, u],
    ])
}

fn vector_cross_prod(a: &Array1<f64>, b: &Array1<f64>) -> Array1<f64> {
    arr1(&[
        -a[[2]] * b[[1]] + a[[1]] * b[[2]],
        -a[[0]] * b[[2]] + a[[2]] * b[[0]],
        -a[[1]] * b[[0]] + a[[0]] * b[[1]],
    ])
}

fn scattering_phase_fn_hg(g: f64, cosq: f64) -> f64 {
    let pfun = 1.0 + g.powi(2) - 2.0 * g * cosq;
    (1.0 - g.powi(2)) / (pfun * pfun.sqrt())
}

// calculate band-mean spectral radiances
fn calc_radiances(config: &Config, dir0: &Array1<f64>, dir1: &Array1<f64>, sinqsol: f64) -> Array1<f64> {
    let waves = &config.waves;
    let mut res = Array1::zeros((waves.len()).f());

    for (iwav, wave) in waves.iter().enumerate() {
        let w = (wave.min + wave.max) / 2.0;
        let rtau = config.input.psfc / 1013.25 * 0.00864 * w.powf(-(3.916 + 0.074 * w + 0.05 / w)); // Rayleigh tau
        let atau = config.input.atau0 * (w / 0.55).powf(-config.input.angexp); // aerosol tau
        let tau_total = rtau + atau;
        let omg = (rtau + config.input.aomg * atau) / tau_total; // single-scattering albedo
        let rsca = rtau / (rtau + config.input.aomg * atau); // Rayleigh scattering fraction
        let asca = config.input.aomg * atau / (rtau + config.input.aomg * atau); // aerosol scattering fraction
        let tau  = tau_total * (1.0 - config.input.chgt);

        // spectral radiance by the single-scattering approximation
        let mu0 = dir0[[2]].abs();
        let cosq = (dir0 * dir1).sum();
        let psi = (1.0 - asca) * 0.75 * (1.0 + cosq.powi(2)) + asca * scattering_phase_fn_hg(config.input.aasy, cosq);

        let (taus, f) = if dir1[[2]] >= 0.0 {
            // upward radiance
            let taus = -(tau_total - tau) * (1.0 / dir1[[2]] + 1.0 / mu0);
            let f = (1.0 - taus.exp()) / (dir1[[2]] + mu0);
            (taus, f)
        } else {
            // downward radiance
            let taus = tau * (1.0 / dir1[[2]] + 1.0 / mu0);
            let f = if taus.abs() < 0.001 {
                // singular case
                -tau / mu0 / dir1[[2]]
            } else {
                // normal case
                (1.0 - taus.exp()) / (dir1[[2]] + mu0)
            };
            (taus, f)
        };

        let acmp = 0.25 * omg * psi * f * (-tau / mu0).exp(); // atmospheric component

        // surface component
        let scmp = if dir1[[2]] >= 0.0 {
            let taus = tau_total / mu0 + (tau_total - tau) / dir1[[2]];
            config.input.galb * taus.exp()
        } else {
            let vp = vector_cross_prod(dir0, dir1);
            let sinq = (&vp * &vp).sum().sqrt();
            if sinq <= sinqsol && cosq > 0.0 {
                let f = 2.0 * PI * sinqsol.powi(2) / (1.0 + (1.0 - sinqsol.powi(2)).sqrt());
                PI / dir0[[2]].abs() / f * (-tau).exp() / dir1[[2]].abs()
            } else {
                0.0
            }
        };

        res.slice_mut(s![iwav]).fill((acmp + scmp) * wave.fsol * dir0[[2]].abs() / PI);
    }

    res
}

fn calc_rimg(config: &Config) -> Array3<f64> {
    let nwav = config.waves.len();
    let mut rimg = Array3::zeros((config.input.nv, config.input.nu, nwav).f());

    // initialize the solar parameters
    let dir0 = new_with_cartesian(-1.0, config.input.the0 * DTOR, config.input.phi0 * DTOR);
    let sinqsol = (config.input.qsol * DTOR).sin();

    // initialize the camera
    let cosq = (config.input.caov.min(179.99) * DTOR).cos();
    let sinq = (config.input.caov.min(179.99) * DTOR).sin();
    let s = 2.0 * sinq * (0.5 / (1.0 + cosq)).sqrt();
    let upmax = s * config.input.nu as f64 / ((config.input.nu as f64).powi(2) + (config.input.cpar * config.input.nv as f64).powi(2)).sqrt();
    let vpmax = s * config.input.cpar * config.input.nv as f64 / ((config.input.nu as f64).powi(2) + (config.input.cpar * config.input.nv as f64).powi(2)).sqrt();
    let rmat = rotate_matrix_z(config.input.cphi * DTOR).dot(&rotate_matrix_y(config.input.cthe * DTOR).dot(&rotate_matrix_z(config.input.cpsi * DTOR)));

    // loop for pixels
    for (iv, mut vv) in rimg.outer_iter_mut().enumerate() {
        for (iu, mut vu) in vv.outer_iter_mut().enumerate() {
            let up = upmax * ((2.0 / config.input.nu as f64) * (iu as f64 + 0.5) - 1.0);
            let vp = vpmax * ((2.0 / config.input.nv as f64) * (iv as f64 + 0.5) - 1.0);
            let s = (up.powi(2) + vp.powi(2)).sqrt();
            let cosq = 1.0 - (s.powi(2) / 2.0);
            let sinq = s / 2.0 * (4.0 - s.powi(2)).sqrt();
            let dir1 = rmat.dot(&arr1(&[
                sinq * up / s,
                sinq * vp / s,
                cosq,
            ]));

            let res = calc_radiances(&config, &dir0, &(-dir1), sinqsol);

            vu.assign(&res);
        }
    }

    rimg
}

fn calc_timg(nv: usize, nu: usize, cimax: f64, waves: &Vec<Wave>, rimg: &Array3<f64>) -> Array3<f64> {
    let mut timg = Array3::zeros((nv, nu, 3).f());
    let nwav = waves.len();

    for (iv, mut vv) in timg.outer_iter_mut().enumerate() {
        for (iu, mut vu) in vv.outer_iter_mut().enumerate() {
            let mut res = Array1::zeros((3).f());

            for (iwav, wave) in waves.iter().enumerate() {
                let cmf = arr1(&[wave.cmf.x, wave.cmf.y, wave.cmf.z]);
                let band_width = wave.max - wave.min;

                res = res + rimg[[iv, iu, iwav]] * cmf * band_width;
            }

            vu.assign(&res)
        }
    }

    let cmat = arr2(&[
                    [3.240970, -1.537383, -0.498611],
                    [-0.969244,  1.875968,  0.041555],
                    [0.055630, -0.203970,  1.056972],
    ]);

    for (iv, mut vv) in timg.outer_iter_mut().enumerate() {
        for (iu, mut vu) in vv.outer_iter_mut().enumerate() {
            let vs = cmat.dot(&vu);
            for (ib, mut vb) in vu.outer_iter_mut().enumerate() {
                let v = (vs[[ib]] / cimax).max(0.0).min(1.0);

                // inverse gamma correction
                if v <= R0 {
                    vb.assign(&arr0(((1.0 + B0) * R0.powf(PWR) - B0) / R0 * v));
                } else {
                    vb.assign(&arr0((1.0 + B0) * v.powf(PWR) - B0));
                }
            }
        }
    }

    timg
}

fn attenuate_if<T: Float + FromPrimitive>(channel: T) -> u8 {
    if channel > FromPrimitive::from_u8(255).unwrap() {
        255u8
    } else if channel < FromPrimitive::from_u8(0).unwrap() {
        0u8
    } else {
        channel.to_u8().unwrap()
    }
}

fn timg_as_u8(timg: &[f64]) -> Vec<u8> {
    timg.iter().map(|x| attenuate_if(*x * 255.99)).collect()
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let config_path: String = args[1].trim().parse().expect("Please provide config.yaml path");

    let config = parse_yaml(&config_path).unwrap();

    validate(&config).unwrap();

    let rimg = calc_rimg(&config);

    let timg = calc_timg(config.input.nv, config.input.nu, config.input.cimax, &config.waves, &rimg);

    let mut timgv = Vec::new();
    for vv in timg.outer_iter() {
        for vu in vv.outer_iter() {
            timgv.push(vu[[0]]);
            timgv.push(vu[[1]]);
            timgv.push(vu[[2]]);
        }
    }

    let timg = timg_as_u8(&timgv);

    image::save_buffer(
        "image.png",
        &timg,
        u32::from_usize(config.input.nu).unwrap(),
        u32::from_usize(config.input.nv).unwrap(),
        image::ColorType::Rgb8,
        ).unwrap()
}
