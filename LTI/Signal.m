classdef Signal < handle
  properties
    m_a;
    m_omega;
    m_lam = 1.0;
    m_phi = 0.0;    
    m_s;
    m_func
  end

  methods
      function this = Signal(options)      
      arguments
          options.r = 1.0; % varying amplitude
          options.omega = 1.0; % frequency
          options.amplitude = 1.0; % initial amplitude
          options.phase = 0.0; % initial phase angle
      end
      this.m_a = log(options.r);
      this.m_omega = options.omega;      
      this.m_lam = options.amplitude * exp(complex(0, options.phase));
      this.m_s = complex(this.m_a, this.m_omega);
      this.m_func = @(t) this.m_lam*exp(this.m_s * t);
    end

    function Show(this)
      n = 1000;
      T = 10;
      t = linspace(0,T,n);
      y = exp(complex(this.m_a, this.m_omega)*t);
      plot(real(y),imag(y));
    end
  end
end