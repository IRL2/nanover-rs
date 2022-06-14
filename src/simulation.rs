extern crate openmm_sys;
use std::cell::RefCell;


use openmm_sys::{
    OpenMM_Force, OpenMM_Force_destroy,
    OpenMM_System,
    OpenMM_System_create, OpenMM_System_destroy,
    OpenMM_System_addForce,
    OpenMM_System_addParticle,
    OpenMM_NonbondedForce,
    OpenMM_NonbondedForce_create, OpenMM_NonbondedForce_destroy,
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
    OpenMM_Context_getPlatform,
    OpenMM_Context_setPositions,
    OpenMM_Context_getState, OpenMM_Context_destroy,
    OpenMM_Platform_getName,
    OpenMM_State_destroy,
    OpenMM_State_DataType_OpenMM_State_Positions,
    OpenMM_State_getTime,
};

pub trait Simulation {
    fn step(&mut self, steps: i32);
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
        println!("Frop ok");
    }
}

impl Simulation for TestSimulation {
    fn step(&mut self, steps: i32) {
        unsafe {println!("{:?}", self.integrator)};
        println!("werp");
        unsafe {
            OpenMM_Integrator_step(self.integrator, steps);
        }
    }
}