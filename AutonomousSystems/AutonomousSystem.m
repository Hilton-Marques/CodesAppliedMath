classdef AutonomousSystem < Scene
properties
	m_F
	m_H
	m_L
	m_G
	m_A
	m_e
	m_x
	m_x_true
end

methods
	function this = AutonomousSystem(F, H , L, G, x0, e0)
		this = this@Scene();
		this.m_F = F;
		this.m_H = H;
		this.m_L = L;
		this.m_G = G;
		this.m_A = (F - L*H);
		this.m_x_true = x0;
		this.m_e = e0;
		this.m_x = x0 - e0;
	end

	function PropagateError(this)
		this.m_e = this.m_A * this.m_e;
	end

	function PropagateTrueSystem(this, t)
		this.m_x_true = this.m_F * this.m_x_true + this.m_G * AutonomousSystem.SamplingInputForce(t);
	end

	function Solver(this)
		delta_t = 0.1;
		n = 50;
		% show initial state
		plot(this.m_x_true(1), this.m_x_true(2), '*','color',this.m_red,'linewidth',1);
		plot(this.m_x(1), this.m_x(2), '*', 'color', this.m_blue,'linewidth',1);
		
		for i = 1:n
			this.PropagateTrueSystem(i*delta_t);
			this.PropagateError();
			this.m_x = this.m_x_true - this.m_e;
			plot(this.m_x_true(1), this.m_x_true(2), '*','color',this.m_red,'linewidth',1);
			plot(this.m_x(1), this.m_x(2), '*', 'color', this.m_blue,'linewidth',1);
			%disp(this.m_e)
			disp(this.m_x_true)
		end
	end
end

methods (Static)
	function f = SamplingInputForce(t)
		f = 1*cos(2*pi*t);
	end
end

end