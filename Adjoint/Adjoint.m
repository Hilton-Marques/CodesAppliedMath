classdef Adjoint < Scene
	properties
		% Adjoint scene properties
	end

	methods
		function this = Adjoint()
			this = this@Scene("adjoint","SE2");
		end

		function ConjugateMap(this, u, Y, X)
			% Conjugate the vector u defined in the tangent space of X
			% to Y. When X is not defined it will be took as the identity. 
			% The formula used is 
			
			% Inputs:
			% 			u: vector to be conjugated
			% 			Y: destiny space
			% 			X: origin space

			arguments
				this
				u (:,1) double
				Y  SE2
				X  SE2 = SE2()
			end

			% Compute the v direction
			v = Y - X;

			% Compute the geodesic from X to Y
			pts = X.geodesic(v, X.m_data);
			n = size(pts,1);

			%Transform u to lie algebra
			U =  X.wedge(u);
			
			% Plot geodeisc
			color = this.m_blue;
			o = X.m_data(1:2,3);
			o = [o;0.1];
			v(3) = o(3);
			this.arrow(o, o + 0.4*v,'tipWidth', 0.015,'color',this.m_colors.MAROON_D,'stemWidth', 0.005);
			coords = zeros(n, 2);
			for i = 1:n	
				color_i = this.m_blue + i/n * (this.m_red - this.m_blue);
      	P_i = squeeze(pts(i,:,:));
				this.drawRobot(P_i, color=color_i);
				o = P_i(1:2,3);
				coords(i,:) = o;
				o = [o;0.1];
				Ux = P_i * U * inv(P_i);
				ui = X.hat(Ux);
				uj = X.hat(P_i * U)
				this.arrow(o, o + ui','tipWidth', 0.015,'color',this.m_colors.GREEN_D,'stemWidth', 0.005);
				this.arrow(o, o + uj','tipWidth', 0.015,'color',this.m_colors.YELLOW_D,'stemWidth', 0.005);
				%quiver(o(1),o(2),ui(1),ui(2));
			end
			this.DrawBoard(n=8);			
			this.setBB(coords');
			camlight
			this.exportFrame(filename="adjoint_se2");
		end

		function LieBracketApprox(this,u, Y, X)
			% Approximate the lie bracket with the adjoint
			% Let du be  difference in velocity in the adjoint path. We have that
			% du = dt * Bracket(U,V), where U is the lie algebra of U and V the
			% lie algebra of Y in X.
			
			% Inputs:
			% 			u: vector to be conjugated
			% 			Y: destiny space
			% 			X: origin space

			arguments
				this
				u (:,1) double
				Y  SE2
				X  SE2 = SE2()
			end

			u1 = u;
			v0 = Y - X;
			V = X.wedge(v0);
			eps = 0.05;
			v = eps*(Y - X);
			%V = X.wedge(v);
			v1 = v;

			X1 = X + v;
			U = X.wedge(u);

			Ux = X1.m_data * U * inv(X1.m_data);
			u2 = X1.hat(Ux)';

			o = X.m_data(1:2,3);
			o = [o;0.1];
			%p = o + v1;
			p = X1.m_data(1:2,3);
			p = [p;0.1];
			q = o + u1;
			r = q + v1;
			v2 = p + u2 - q;
			bracket = X1.hat(V*U - U*V)';
			du = eps * bracket;
			v3 = v1 + du;
			
			% Main approximation
			brack_approx = (u2 - u1)/eps;

			fprintf('the error is in the bracket approximation is %d \n', norm(bracket - brack_approx));

			this.arrow(o, o + u1,'tipWidth', 0.010,'color',this.m_colors.GREEN_B,'stemWidth', 0.005);
			this.arrow(o, o + v1,'tipWidth', 0.010,'color',this.m_colors.MAROON_B,'stemWidth', 0.005);
			this.arrow(p, p + u2,'tipWidth', 0.010,'color',this.m_colors.GREEN_C,'stemWidth', 0.005);
			this.arrow(q, q + v2,'tipWidth', 0.010,'color',this.m_colors.RED_D,'stemWidth', 0.005);
			this.arrow(q, q + v1,'tipWidth', 0.010,'color',this.m_colors.MAROON_B,'stemWidth', 0.005);
			this.arrow(r, r + du,'tipWidth', 0.005,'color',this.m_colors.YELLOW_D,'stemWidth', 0.002);

% 			quiver(o(1), o(2), u1(1), u1(2),'off');
% 			quiver(o(1), o(2), v1(1), v1(2),'off');			
% 			quiver(p(1), p(2), u2(1), u2(2),'off');
% 			quiver(q(1), q(2), v2(1), v2(2),'off');
% 			quiver(q(1), q(2), v1(1), v1(2),'off');
% 			quiver(r(1), r(2), du(1), du(2),'off');

			this.drawRobot(X.m_data, color=this.m_blue, s = 0.15);
			
			color = this.m_blue + eps*(this.m_red - this.m_blue);
			this.drawRobot(X1.m_data, color=color, s = 0.15);
			o = [0;0;0.1];
			v0(3) = 0.1;

			%this.arrow(o, o + 0.4*v0,'tipWidth', 0.015,'color',this.m_colors.MAROON_D,'stemWidth', 0.005);
			this.setBB([[0,0,0];[0.8, 0.4,0]]',0.2);
			this.DrawBoard(n=8);
			camlight
			this.exportFrame(filename="bracket_se2");			
		end

		function FrameTransformation(this, A, B)
			I = SE2();
			this.drawRobot(I.m_data, color=this.m_green);
			this.drawRobot(A.m_data, color=this.m_blue);
			this.drawRobot(B.m_data, color=this.m_red);

			C = B ^ A;
			D = A.inverse ^ B ^ A;
			v = D.log();
			u = B.log();
			w = A.log();
			U = A.wedge(u);
			W = inv(A.m_data) * A.wedge(u) * A.m_data;
			v_ = A.hat(W)';
			assert(norm(v_ - v) < 1e-8);
			%E = A ^ B;
% 			check = B.m_data(1:2,1:2)*A.m_data(1:2,3) + B.m_data(1:2,3) - A.m_data(1:2,3);
% 			check2 = D.m_data(1:2, 3);
% 			norm(check - check2)

			this.drawRobot(C.m_data, color=this.m_colors.MAROON_A);
			this.drawRobot(D.m_data, color=this.m_colors.MAROON_D);

			o = [0;0;0.0];
			fac = 0.6;

			this.setBB([[-0.8,0,0];...
									[C.m_data(1,3), C.m_data(2,3),0]; ...
									[A.m_data(1,3), A.m_data(2,3),0]; ...
									[D.m_data(1,3), D.m_data(2,3),0]]',0.33);
			this.DrawBoard(n=8);

			this.arrow(o, o + fac*v,'tipWidth', 0.035,'color',this.m_colors.MAROON_D,'stemWidth', 0.015);
			this.arrow(o, o + fac*u,'tipWidth', 0.035,'color',this.m_red,'stemWidth', 0.015);
			w(3) = 0.0;
		  this.arrow(o, o + fac*w,'tipWidth', 0.035,'color',this.m_blue,'stemWidth', 0.015);
			camlight;

			%Show geodesics
			pts = A.geodesic(A.log(), eye(3),linspace(0,1,4));

			n = size(pts,1);			
			I = SE2();
			for i = 1:n				
				color_i = this.m_red + i/n * (this.hex2rgb(this.m_colors.MAROON_D) - this.m_red);
      	P_i = squeeze(pts(i,:,:));
				Ux = inv(P_i) * U * P_i;
				u_i = A.hat(Ux)';
				R_i = I + u_i;
				this.drawRobot(R_i.m_data, color=color_i);
				this.arrow(o, o + fac*u_i,'tipWidth', 0.035,'color',color_i,'stemWidth', 0.015);
			end

			B_algebra = B.wedge(B.log());
			A_algebra = A.wedge(A.log());
			bracket = A.hat(B_algebra*A_algebra - A_algebra*B_algebra)';			
			bracket2 = A.hat(A_algebra*W - W*A_algebra)';
			%this.arrow(fac*u, fac*u + fac*0.2*bracket,'tipWidth', 0.02,'color',this.m_colors.YELLOW_D,'stemWidth', 0.012);
			this.arrow(fac*v, fac*v + fac*0.4*bracket2,'tipWidth', 0.02,'color',this.m_colors.YELLOW_D,'stemWidth', 0.012);
			this.arrow(o, o + fac*bracket2,'tipWidth', 0.035,'color',this.m_colors.YELLOW_D,'stemWidth', 0.015);

			%this.drawRobot(E.m_data, color=this.m_colors.PINK);
			this.exportFrame(filename="action_IV");
		end
		function error = adjoint_rep(this, X,Y)
			x = X.log();
			X_algebra = X.wedge(x);
			v = Y.log();
			U = inv(X.m_data) * X.wedge(v) * X.m_data;
			real = v - X.hat(U)';
			approx = X.hat(this.compute_conjugate_difference(X_algebra, U,6))';
			%% Why the following gives the exact value for the above approximation?
			%% This works for all Lie Groups?
			bracket_exact = X.hat(X.m_data*U - U*X.m_data)';
			error = norm(real - approx);
			error = norm(real - bracket_exact);
			keyboard;
		end
	end
	methods (Static)
		function result = compute_conjugate_difference(X, Y, order)
			% Compute U - Y where U = XYX^-1 up to a specified order of approximation
			% Inputs:
			%   X - matrix representing element X in the Lie algebra
			%   Y - matrix representing element Y in the Lie algebra
			%   order - the order of approximation for the series expansion
			% Output:
			%   result - matrix representing the difference U - Y

			% Initialize result with the first term in the series
			result = X * Y - Y * X;

			% Initialize the current term with the first commutator
			current_term = result;

			% Loop to compute higher-order commutators up to the specified order
			for n = 2:order
				current_term = X * current_term - current_term * X; % Compute the next commutator
				result = result + (1/factorial(n)) * current_term; % Add the next term in the series
			end
		end
	end
end