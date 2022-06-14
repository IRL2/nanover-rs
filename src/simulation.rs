extern crate openmm_sys;

use openmm_sys::{
    OpenMM_Force,
    OpenMM_System,
    OpenMM_System_create, OpenMM_System_destroy,
    OpenMM_System_addForce,
    OpenMM_System_addParticle,
    OpenMM_NonbondedForce_create,
    OpenMM_NonbondedForce_addParticle,
    OpenMM_Vec3,
    OpenMM_Vec3Array,
    OpenMM_Vec3Array_create, OpenMM_Vec3Array_destroy,
    OpenMM_Vec3Array_set,
    OpenMM_Integrator, OpenMM_Integrator_destroy,
    OpenMM_Integrator_step,
    OpenMM_VerletIntegrator_create,
    OpenMM_Context,
    OpenMM_Context_create,
    OpenMM_Context_setPositions,
    OpenMM_Context_destroy,
    OpenMM_Context_getState,
    OpenMM_State_destroy,
    OpenMM_State_DataType_OpenMM_State_Positions,
    OpenMM_State_getPositions,
    OpenMM_Vec3Array_getSize,
    OpenMM_Vec3Array_get,
    OpenMM_Vec3_scale,
};

use crate::frame::FrameData;

pub trait Simulation {
    fn step(&mut self, steps: i32);
}

pub trait ToFrameData {
    fn to_framedata(&self) -> FrameData;
}

pub struct TestSimulation {
    system: *mut OpenMM_System,
    init_pos: *mut OpenMM_Vec3Array,
    integrator: *mut OpenMM_Integrator,
    context: *mut OpenMM_Context,
}

impl TestSimulation {
    pub fn new() -> Self {
        let sim = unsafe {
            let system = OpenMM_System_create();
            let nonbond = OpenMM_NonbondedForce_create();
            OpenMM_System_addForce(system, nonbond as *mut OpenMM_Force);
            let init_pos = OpenMM_Vec3Array_create(3);

            for a in 0..3 {
                let pos = OpenMM_Vec3 {
                    x: 0.95 * a as f64,
                    y: 0.95 * a as f64,
                    z: 0.95 * a as f64,
                };

                OpenMM_Vec3Array_set(init_pos, a, pos);
                OpenMM_System_addParticle(system, 39.95);
                OpenMM_NonbondedForce_addParticle(nonbond, 0.0, 0.3350, 0.996);
            }

            let integrator = OpenMM_VerletIntegrator_create(0.004) as *mut OpenMM_Integrator;
            let context = OpenMM_Context_create(system, integrator);
            OpenMM_Context_setPositions(context, init_pos);
            Self {system, init_pos, integrator, context}
        };
        sim
    }
}

impl Drop for TestSimulation {
    fn drop(&mut self) {
        println!("Frop!");
        unsafe {
            OpenMM_Vec3Array_destroy(self.init_pos);
            OpenMM_Context_destroy(self.context);
            OpenMM_Integrator_destroy(self.integrator);
            OpenMM_System_destroy(self.system);
        }
    }
}

impl Simulation for TestSimulation {
    fn step(&mut self, steps: i32) {
        unsafe {
            OpenMM_Integrator_step(self.integrator, steps);
        }
    }
}

impl ToFrameData for TestSimulation {
    fn to_framedata(&self) -> FrameData {
        let mut positions = Vec::<f32>::new();
        unsafe {
            let state = OpenMM_Context_getState(
                self.context,
                OpenMM_State_DataType_OpenMM_State_Positions as i32,
                0,
            );
            let pos_state = OpenMM_State_getPositions(state);
            let particle_count = OpenMM_Vec3Array_getSize(pos_state);
            println!("{particle_count}");
            for i in 0..particle_count {
                let pos = OpenMM_Vec3_scale(*OpenMM_Vec3Array_get(pos_state, i), 1.0);
                positions.push(pos.x as f32);
                positions.push(pos.y as f32);
                positions.push(pos.z as f32);
            }
            OpenMM_State_destroy(state);
        }
        let mut frame = FrameData::empty();
        println!("{:?}", positions);
        frame.insert_number_value("particle.count", (positions.len() / 3) as f64).unwrap();
        frame.insert_float_array("particle.positions", positions).unwrap();
        frame
    }
}