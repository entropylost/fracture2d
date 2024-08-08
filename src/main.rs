// Version of the fracture simulation from "Simulating Fractures with Bonded Discrete Element Method" paper
// Taken from the supplementary material and translated to Rust.

use std::f64::consts::PI;

use macroquad::prelude::*;

#[derive(Debug, Clone, Copy)]
struct Bond {
    broken: bool,
    endpoint: u32,
    length: f64,
    direction: DVec2,
    max_normal_force: f64,
    max_tangent_force: f64,
}
impl Bond {
    fn new(j: u32, l: f64, d: DVec2, tn: f64, tt: f64) -> Self {
        Self {
            endpoint: j,
            broken: false,
            length: l,
            direction: d,
            max_normal_force: tn,
            max_tangent_force: tt,
        }
    }
}

struct Particle {
    bonds: Vec<Bond>,
    inverse_mass: f64,
    inverse_moment: f64,
    radius: f64,
    angle: f64,
    angvel: f64,
    // For velocity verlet.
    angvel_mid: f64,
    torque: f64,
    position: DVec2,
    velocity: DVec2,
    // For velocity verlet.
    velocity_mid: DVec2,
    force: DVec2,
    color: Color,
}
impl Particle {
    fn new(im: f64, ii: f64, r: f64, x: DVec2, c: Color) -> Self {
        Self {
            bonds: vec![],
            inverse_mass: im,
            inverse_moment: ii,
            radius: r,
            angle: 0.0,
            angvel: 0.0,
            angvel_mid: 0.0,
            torque: 0.0,
            position: x,
            velocity: DVec2::ZERO,
            velocity_mid: DVec2::ZERO,
            force: DVec2::ZERO,
            color: c,
        }
    }
}

fn range(start: f64, end: f64, step: f64) -> impl Iterator<Item = f64> {
    let mut x = start;
    std::iter::from_fn(move || {
        if x < end {
            let y = x;
            x += step;
            Some(y)
        } else {
            None
        }
    })
}

fn clamp_angle(a: f64) -> f64 {
    (a + PI).rem_euclid(2.0 * PI) - PI
}

#[macroquad::main("Fracture 2d")]
async fn main() {
    let fps: f64 = 60.0;
    let r: f64 = 0.02;
    let r2 = r * r;
    let r4 = r2 * r2;
    // Young's modulus
    let kn: f64 = 1e7;
    let particle_invmass = 1.3e-4 / r2;
    let particle_invmoment = 2.6e-4 / r4;
    let mut pts = vec![];
    for x in range(0.1, 0.9, 2.0 * r) {
        for y in range(0.4, 0.45, 2.0 * r) {
            pts.push(Particle::new(
                particle_invmass,
                particle_invmoment,
                r,
                DVec2::new(x, y),
                GREEN,
            ));
        }
    }
    for x in range(0.7, 0.9, 2.0 * r) {
        for y in range(0.1, 0.3, 2.0 * r) {
            pts.push(Particle::new(
                particle_invmass,
                particle_invmoment,
                r,
                DVec2::new(x, y),
                BLUE,
            ));
        }
    }
    for i in 0..pts.len() {
        for j in 0..pts.len() {
            if i == j {
                continue;
            }
            let l = pts[j].position - pts[i].position;
            let overlap = pts[i].radius + pts[j].radius - l.length();
            if overlap >= -0.1 * r {
                pts[i].bonds.push(Bond::new(
                    j as u32,
                    2.0 * r,
                    l.normalize(),
                    0.07 * kn,
                    0.07 * kn,
                ));
            }
        }
    }
    for x in range(r, 1.0, 2.0 * r) {
        pts.push(Particle::new(0.0, 0.0, r, DVec2::new(x, 0.0), GRAY));
        pts.push(Particle::new(0.0, 0.0, r, DVec2::new(x, 1.0), GRAY));
        pts.push(Particle::new(0.0, 0.0, r, DVec2::new(0.0, x), GRAY));
        pts.push(Particle::new(0.0, 0.0, r, DVec2::new(1.0, x), GRAY));
    }
    for x in range(0.2, 0.8, 2.0 * r) {
        for y in range(0.5, 0.7, 2.0 * r) {
            pts.push(Particle::new(
                particle_invmass,
                particle_invmoment,
                r,
                DVec2::new(x, y),
                ORANGE,
            ));
        }
    }
    // Extra timestep?
    let s = (1.0 / fps / (7.5e3 * r2 / kn)) as u32 * 10;
    let dt = 1.0 / fps / s as f64;
    println!("S: {:?}", s);
    loop {
        for _ in 0..1000 {
            for particle in &mut pts {
                particle.force = DVec2::ZERO;
                particle.torque = 0.0;
                particle.position += particle.velocity_mid * dt;
                particle.angle += particle.angvel_mid * dt;
                particle.angle = clamp_angle(particle.angle);
            }
            for i in 0..pts.len() {
                if pts[i].inverse_mass <= 1e-6 {
                    continue;
                }
                for j in 0..pts.len() {
                    if i == j {
                        continue;
                    }
                    let lij = pts[i].position - pts[j].position;
                    let o = pts[i].radius + pts[j].radius - lij.length();
                    if o < 1e-12 {
                        continue;
                    }
                    let n = lij.normalize();
                    let a = 1.4 * (kn / (pts[i].inverse_mass + pts[j].inverse_mass)).sqrt();
                    let force = n * (kn * o + a * (pts[j].velocity - pts[i].velocity).dot(n));
                    pts[i].force += force;
                }
                for iter in 0..pts[i].bonds.len() {
                    let b = &pts[i].bonds[iter];
                    if b.broken {
                        continue;
                    }
                    let l = pts[b.endpoint as usize].position - pts[i].position;
                    let n = l.normalize();
                    let t = DVec2::new(-n.y, n.x);
                    let dl = l.length() - b.length;
                    let qb = b.direction.y.atan2(b.direction.x) - n.y.atan2(n.x);
                    let ti = clamp_angle(qb + pts[i].angle);
                    let tj = clamp_angle(qb + pts[b.endpoint as usize].angle);
                    let f_n = n * kn * dl;
                    let f_t = t * -kn / 3.0 * r2 / l.length() * (ti + tj);
                    let t = kn / 6.0 * r2 * (tj - 3.0 * ti);
                    if (dl > 0.0
                        && (f_n.length() / 2.0 / r + (kn / 2.0 * (tj - ti)).abs())
                            > b.max_normal_force)
                        || f_t.length() / 2.0 / r > b.max_tangent_force
                    {
                        pts[i].bonds[iter].broken = true;
                        continue;
                    }
                    pts[i].force += f_n + f_t;
                    pts[i].torque += t;
                }
                for particle in &mut pts {
                    let acc = particle.force * particle.inverse_mass
                        + if particle.inverse_mass > 1e-6 {
                            DVec2::new(0.0, -9.8)
                        } else {
                            DVec2::ZERO
                        };
                    particle.velocity = particle.velocity_mid + acc * 0.5 * dt;
                    particle.velocity_mid += acc * dt;
                    particle.angvel =
                        particle.angvel_mid + particle.torque * particle.inverse_moment * 0.5 * dt;
                    particle.angvel_mid += particle.torque * particle.inverse_moment * dt;
                }
            }
        }

        // Rendering

        let scaling = 500.0;

        clear_background(WHITE);
        for p in &pts {
            let a = p.position - DVec2::splat(0.5);

            draw_circle(
                screen_height() / 2.0 + a.x as f32 * scaling,
                screen_width() / 2.0 - a.y as f32 * scaling,
                p.radius as f32 * scaling,
                p.color,
            );

            for b in &p.bonds {
                if !b.broken {
                    continue;
                }
                let b = pts[b.endpoint as usize].position - DVec2::splat(0.5);

                draw_line(
                    screen_height() / 2.0 + a.x as f32 * scaling,
                    screen_width() / 2.0 - a.y as f32 * scaling,
                    screen_height() / 2.0 + b.x as f32 * scaling,
                    screen_width() / 2.0 - b.y as f32 * scaling,
                    3.0,
                    BLACK,
                );
            }
        }
        next_frame().await
    }
}
