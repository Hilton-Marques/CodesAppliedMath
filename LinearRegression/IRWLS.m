classdef IRWLS < handle
	properties
		m_H
		m_b
		m_psi
		m_iter
		m_rho
		m_x0
		m_tol = 1e-5
		m_M
	end

	methods
		function this = IRWLS(H, b, method, M, x0, iters)
			arguments
				H (:,:) double
				b (:,1) double
				method char
				M = []
				x0 = []
				iters (1,1) double = 10				
			end
			%this = this@Scene();
			if (isempty(x0))
				x0 = pinv(H)*b;
			end
			this.m_x0 = x0;
			this.m_H = H;
			this.m_b = b;
			this.m_iter = iters;
			switch method
				case 'l2'
					this.m_psi = @(x) ones(size(x));
					this.m_rho = @(x) x.^2;
				case 'huber'
					%this value is defined in zoubir , pg.23
					are_huber = M;
					if (isempty(are_huber))
						are_huber = 1.345;
					end
					this.m_rho = @(x) this.huberRho(x, are_huber);
					this.m_psi = @(x) this.huberPsi(x, are_huber);
				case 'tukey'
					are_tukey = M;
					if (isempty(are_tukey))
						are_tukey = 3.4437;
					end
					this.m_rho = @(x) this.tukeyRho(x, are_tukey);
					this.m_psi = @(x) this.tukeyPsi(x, are_tukey);
				otherwise
					error('Unknown method');
			end
		end

		function x = Solve(this)
			% Solve the problem
			x = this.m_x0;
			% Compute the residuals
			for i = 1:this.m_iter
				r = abs(this.m_b - this.m_H*x);
				% Avoid division by zero
				r(r<this.m_tol) = this.m_tol;
				% Compute the weights
				W = this.m_psi(r);
				% Compute the weighted least squares solution
				S = this.m_H' * diag(W);
				x_n = S*this.m_H \ S*this.m_b;
				e = x_n - x;
				er = dot(e,e)/dot(x,x);
				if (sqrt(er) < this.m_tol)
					x = x_n;
					return;
				end
				x = x_n;
			end
		end

		function z = Show(this, x)
			% Create original function
			r = this.m_H*x - this.m_b;
			z_value = this.J(r);

			% Create Contour plot
			n = 300;
			x_inv = linspace(-10,10,n);
			y_inv = linspace(-10,10,n);
			[x_grid, y_grid] = meshgrid(x_inv,y_inv);
			x = [x_grid(:), y_grid(:)];
			values = this.J(x);
			z = reshape(values, n, n);
			[z,h] = contourf(x_grid, y_grid, z, [z_value,z_value],'FaceAlpha',0.7);
			delete(h);
			z = z(:,2:end);
		end
		function res = J(this, r)
			% Compute the objective function
			if size(r,2) == 1
				r = r';
			end
			n = size(r,1);
			res = zeros(n,1);
			for i = 1:n
				res(i) = sum(this.m_rho(r(i,:)));
			end
		end
	end

	methods (Static)
		function x = huberRho(x, c)
			% Huber cost
			x = abs(x);
			x(x <= c) = 0.5*x(x <= c).^2;
			x(x > c) = c*x(x > c) - 0.5*c^2;
		end

		function x = huberPsi(x, c)
			% Huber Psi (Zoubir, pg. 12)
			x_abs = abs(x);
			x(x_abs <= c) = 1;
			x(x_abs > c) = c ./ x_abs((x_abs > c));
		end

		function x = tukeyRho(x, c)
			% Tukey function
			x = abs(x);
			x(x <= c) = (c^2/6)*(1 - (1 - (x(x <= c)/c).^2).^3);
			x(x > c) = c^2/6;
		end

		function x = tukeyPsi(x, c)
			% Tukey Psi (Zoubir, pg. 12)
			x = abs(x);
			x(x <= c) = (1 - (x(x <= c)/c).^2).^2;
			x(x > c) = 0;
		end
	end
end