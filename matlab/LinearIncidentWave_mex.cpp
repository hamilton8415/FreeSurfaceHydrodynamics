#include <iostream>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "FreeSurfaceHydrodynamics/LinearIncidentWave.hpp"

class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) override {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;

        if (inputs.empty() || inputs[0].getType() != matlab::data::ArrayType::CHAR) {
            //matlabPtr->error("First argument must be an action string ('New', 'Delete').");
            return;
        }

        // Convert the action argument to a standard C++ string
        matlab::data::CharArray actionArr = inputs[0];
        std::string action = actionArr.toAscii();

        // ----------------- ACTION: CREATE -----------------
        if (action == "New") {


            //int initVal = (inputs.size() > 1) ? inputs[1][0] : 0.0;
            
            // Instantiate the C++ class on the heap
            switch(inputs.size())
            {
              case 1:  //No additional arguments
              {
                LinearIncidentWave* obj = new LinearIncidentWave();
                // Cast the pointer to a uint64 to pass safely to MATLAB
                matlab::data::TypedArray<uint64_t> handle = factory.createScalar<uint64_t>(reinterpret_cast<uint64_t>(obj));
                outputs[0] = handle;
                break;
              }
              case 2:  // One additional argument
              {
                unsigned int newVal = inputs[2][0];
                LinearIncidentWave* obj = new LinearIncidentWave(newVal);
                // Cast the pointer to a uint64 to pass safely to MATLAB
                matlab::data::TypedArray<uint64_t> handle = factory.createScalar<uint64_t>(reinterpret_cast<uint64_t>(obj));
                outputs[0] = handle;
                break;
              }
            }
        }

           
        // ----------------- SetSeed(unsigned int seed) -----------------
        else if (action == "SetSeed") {
            uint64_t handle = inputs[1][0];
            unsigned int newVal = inputs[2][0];
            
            LinearIncidentWave* obj = reinterpret_cast<LinearIncidentWave*>(handle);
            obj->SetSeed(newVal);
        }

        // ----------------- SetToMonoChromatic -----------------
       else if (action == "SetToMonoChromatic") {
            uint64_t handle = inputs[1][0];
            double A = inputs[2][0];
            double T = inputs[3][0];
            double phase = inputs[4][0];
            double beta = inputs[5][0];
            
            LinearIncidentWave* obj = reinterpret_cast<LinearIncidentWave*>(handle);
            obj->SetToMonoChromatic(A,T,phase,beta);

       }

      
        // ----------------- eta -----------------
       else if (action == "eta") {
          uint64_t handle = inputs[1][0];
          LinearIncidentWave* obj = reinterpret_cast<LinearIncidentWave*>(handle);
          switch(inputs.size())  // First argument is name, second argument is handle, then the function arguments
          {
            case 5:  // Three arguments  eta(double x, double y, double t)
              {
               double x = inputs[2][0];
               double y = inputs[3][0];
               double t = inputs[4][0];
               switch(outputs.size())
               {
                 case 1:
                   outputs[0] = factory.createScalar<double>(obj->eta(x,y,t));
                 break;

                 case 3:
                 double deta_dx;
                 double deta_dy;
                 double val = obj->eta(x,y,t,&deta_dx,&deta_dy); 
                   outputs[0] = factory.createScalar<double>(val);
                   outputs[1] = factory.createScalar<double>(deta_dx);
                   outputs[2] = factory.createScalar<double>(deta_dy);
                 break;
               }
              }
              break;
            case 6:  // Four arguments  eta(double x, double y, double t, int n
              { 
              double x = inputs[2][0];
              double y = inputs[3][0];
              double t = inputs[4][0];
              int n = inputs[5][0];
              outputs[0] = factory.createScalar<double>(obj->eta(x,y,t,n));
              }
              break;
          }
            
       }
       


        
        // ----------------- Int MajorVersionNumber(); -----------------
        else if (action == "MajorVersionNumber") {
            uint64_t handle = inputs[1][0];
            LinearIncidentWave* obj = reinterpret_cast<LinearIncidentWave*>(handle);
            
            int val = obj->MajorVersionNumber();
            outputs[0] = factory.createScalar<int>(val);
        }

        // ----------------- Int MinorVersionNumber(); -----------------
        else if (action == "MinorVersionNumber") {
            uint64_t handle = inputs[1][0];
            LinearIncidentWave* obj = reinterpret_cast<LinearIncidentWave*>(handle);
            
            int val = obj->MinorVersionNumber();
            outputs[0] = factory.createScalar<int>(val);
        }
        
        // ----------------- Int PatchVersionNumber(); -----------------
        else if (action == "PatchVersionNumber") {
            uint64_t handle = inputs[1][0];
            LinearIncidentWave* obj = reinterpret_cast<LinearIncidentWave*>(handle);
            
            int val = obj->PatchVersionNumber();
            outputs[0] = factory.createScalar<int>(val);
        }
        
        // ----------------- ACTION: DELETE -----------------
        else if (action == "Delete") {
            if (inputs.size() < 2) return;
            uint64_t handle = inputs[1][0];
            LinearIncidentWave* obj = reinterpret_cast<LinearIncidentWave*>(handle);
            
            // Clean up heap memory safely
            delete obj;
        }
        
        else {
            //matlabPtr->error("Unknown action string provided.");
        }
    }
};