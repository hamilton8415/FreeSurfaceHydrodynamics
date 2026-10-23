classdef LinearIncidentWave < handle
    properties (Access = private)
        ObjectHandle % Holds the uint64 representation of the C++ pointer
    end
    
    methods
        % Constructor
        function obj = LinearIncidentWave()
            % Call MEX to instantiate C++ object
            obj.ObjectHandle = LinearIncidentWave_mex('New');
        end

        %SetSeed method
        function SetSeed(obj, val)
            LinearIncidentWave_mex('SetSeed', obj.ObjectHandle, val);
        end

        %SetToMonoChromatic method
        function SetToMonoChromatic(obj, A, T, phase, beta)
            LinearIncidentWave_mex('SetToMonoChromatic', obj.ObjectHandle, A, T, phase, beta);
        end


        %SetToPiersonMoskowitzSpectrum method
        function SetToPiersonMoskowitzSpectrum(obj, varargin)
            if(nargin == 3)    % (double Hs, double beta)
              LinearIncidentWave_mex('SetToPiersonMoskowitzSpectrum', obj.ObjectHandle, ...
                                     varargin{1}, varargin{2});
            end
            if(nargin == 4)    % (double Hs, double beta, int n_phases)
              LinearIncidentWave_mex('SetToPiersonMoskowitzSpectrum', obj.ObjectHandle, ...
                                     varargin{1}, varargin{2}, varargin{3});
             end
         end

        %SetToPiersonMoskowitzSpectrumWithCos2Spreading method
        function SetToPiersonMoskowitzSpectrumWithCos2Spreading(obj, varargin)
            if(nargin == 5)    % (double Hs, double beta,  int spreading_factor, int n_sectors)
              LinearIncidentWave_mex('SetToPiersonMoskowitzSpectrumWithCos2Spreading', obj.ObjectHandle, ...
                                      varargin{1}, varargin{2}, varargin{3}, varargin{4});
            end
            if(nargin == 6)    % (double Hs, double beta, int n_phases  int spreading_factor, int n_sectors)
              LinearIncidentWave_mex('SetToPiersonMoskowitzSpectrumWithCos2Spreading', obj.ObjectHandle, ...
                                      varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
             end
         end

        %SetToBretschneiderSpectrum method
        function SetToBretschneiderSpectrum(obj, varargin)
            if(nargin == 4)    % (double Hs, double beta)
              LinearIncidentWave_mex('SetToBretschneiderSpectrum', obj.ObjectHandle, ...
                                     varargin{1}, varargin{2}, varargin{3});
            end
            if(nargin == 5)    % (double Hs, double beta, int n_phases)
              LinearIncidentWave_mex('SetToBretschneiderSpectrum', obj.ObjectHandle, ...
                                     varargin{1}, varargin{2}, varargin{3}, varargin{4});
             end
         end

        %SetToBretschneiderSpectrumWithCos2Spreading method
        function SetToPiersonMoskowitzSpectrumWithCos2Spreading(obj, varargin)
            if(nargin == 6)    % (double Hs, double beta,  int spreading_factor, int n_sectors)
              LinearIncidentWave_mex('SetToBretschneiderSpectrumWithCos2Spreading', obj.ObjectHandle, ...
                                      varargin{1}, varargin{2}, varargin{3}, varargin{4}), varargin{5},;
            end
            if(nargin == 7)    % (double Hs, double beta, int n_phases  int spreading_factor, int n_sectors)
              LinearIncidentWave_mex('SetToBretschneiderSpectrumWithCos2Spreading', obj.ObjectHandle, ...
                                      varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5}, varargin{6},);
             end

        % eta method
        function varargout = eta(obj,varargin)
             if(nargout == 1)
               if(nargin == 4)
                 varargout{1} = LinearIncidentWave_mex('eta', obj.ObjectHandle, varargin{1}, varargin{2}, varargin{3});
               end
               if(nargin == 5)
                 varargout{1} = LinearIncidentWave_mex('eta', obj.ObjectHandle, varargin{1}, varargin{2}, varargin{3}, varargin{4});
               end
             end
             if(nargout == 3)
               if(nargin == 4)
                 [varargout{1} varargout{2} varargout{3}] = LinearIncidentWave_mex('eta', obj.ObjectHandle, varargin{1}, varargin{2}, varargin{3});
               end
               if(nargin == 5)
                 [varargout{1}  varargout{2} varargout{3}] = LinearIncidentWave_mex('eta', obj.ObjectHandle, varargin{1}, varargin{2}, varargin{3}, varargin{4});
               end
             end
             if(nargout == 6)
               if(nargin == 4)
                 [varargout{1} varargout{2} varargout{3} varargout{4} varargout{5}] = LinearIncidentWave_mex('eta', obj.ObjectHandle, varargin{1}, varargin{2}, varargin{3});
               end
               if(nargin == 5)
                 [varargout{1}  varargout{2} varargout{3} varargout{4} varargout{5}] = LinearIncidentWave_mex('eta', obj.ObjectHandle, varargin{1}, varargin{2}, varargin{3}, varargin{4});
               end
             end
        end

        % etadot method
        function varargout = etadot(obj,varargin)
             %if(nargout == 1)
               if(nargin == 4)
                 [varargout{1}] = LinearIncidentWave_mex('etadot', obj.ObjectHandle, varargin{1}, varargin{2}, varargin{3});
               end
               if(nargin == 5)
                 [varargout{1}] = LinearIncidentWave_mex('etadot', obj.ObjectHandle, varargin{1}, varargin{2}, varargin{3}, varargin{4});
               end
             %end
        end

        % MajorVersionNumber method
        function val = MajorVersionNumber(obj)
            val = LinearIncidentWave_mex('MajorVersionNumber', obj.ObjectHandle);
        end

        % MinorVersionNumber method
        function val = MinorVersionNumber(obj)
            val = LinearIncidentWave_mex('MinorVersionNumber', obj.ObjectHandle);
        end

        % PatchVersionNumber method
        function val = PatchVersionNumber(obj)
            val = LinearIncidentWave_mex('PatchVersionNumber', obj.ObjectHandle);
        end 
        
        % Destructor (triggers automatically when MATLAB clears the object)
        function delete(obj)
            if ~isempty(obj.ObjectHandle)
                LinearIncidentWave_mex('Delete', obj.ObjectHandle);
                obj.ObjectHandle = [];
            end
        end
    end
end
