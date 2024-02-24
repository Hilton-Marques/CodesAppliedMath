classdef QRGivensRotations < Scene
    properties
        m_A
        m_handles = [];
    end
    methods
        function this = QRGivensRotations(A)
            this.m_A = A;
            this.setBB([A,[0; 0; 0]]);
            this.show();
            %%% Scene 1 - Proj zy %%%
            this.reset();
            v = this.proj([1;0;0], A(:,1), 100);
            this.m_handles = [this.m_handles, arrow3([0,0,0], v', 'v')];
            %this.save('projection.gif');

            %%% Sceen 2 - Rotation every by Givens rotations %%%%
            this.reset();
            G = eye(3);
            [c,s] = this.givensrotation(A(2,1), A(3,1));
            G([2, 3],[2, 3]) = [c -s; s c];
            this.rot(G', [A,v])
            %this.save('first_rotation.gif');
            A = G' * A;
            %%% Scene 3 - Full Givens %%%%
            this.qrgivens(A);
        end
        function [Q,R] = qrgivens(this, A)
            [m,n] = size(A);
            Q = eye(m);
            R = A;
            for j = 1:n
                for i = m:-1:(j+1)
                    G = eye(m);
                    this.reset();
                    [c,s] = this.givensrotation( R(i-1,j),R(i,j) );
                    if (c == 1)
                        continue;
                    end
                    G([i-1, i],[i-1, i]) = [c -s; s c];
                    %this.rot(G', R);
                    %this.save(strcat(num2str(i),'_rotation.gif'));
                    R = G'*R;                   
                    Q = Q*G;
                end
            end
        end
        %proj v in plane given by n
        function vi = proj(this, n, v, n_frames)
            time = 1;
            n = n/norm(n);
            A = (eye(3) - n * n');
            for i = 1:n_frames
                I = eye(3) + (i/n_frames)*(A - eye(3));
                vi = I*v;
                h = arrow3([0,0,0], vi', 'v');
                this.get();
                delete(h);
            end
            vi = A*v;
        end
        %rotate by R matrix a set of P points
        function r = rot(this, R, P)
            colors = ['k', 'k', 'v'];
            [m, n] = size(P);
            eul = rotm2eul(R);
            n_frames = 100;
            for i = 1:n_frames
                delete(this.m_handles);
                InAngle = rotm2eul(eye(3)) + (i/n_frames)*(eul-zeros(1,3));
                InRot = eul2rotm(InAngle);
                IN_P = InRot * P;
                h = [];
                for j = 1:n
                    h = [h, arrow3([0,0,0], IN_P(:,j)', colors(j))];
                end
                this.setBB([IN_P, eye(3)]);
                this.get();
                this.m_handles = h;
            end
        end
        function [c,s] = givensrotation(this, a,b)
            if b == 0
                c = 1;
                s = 0;
            else
                if abs(b) > abs(a)
                    r = a / b;
                    s = 1 / sqrt(1 + r^2);
                    c = s*r;
                else
                    r = b / a;
                    c = 1 / sqrt(1 + r^2);
                    s = c*r;
                end
            end
        end
        function show(this)
            view(126, 26);
            this.m_handles = arrow3([0 0 0], this.m_A(:,1)');
            this.m_handles = [this.m_handles, arrow3([0 0 0], this.m_A(:,2)')];
            arrow3([0 0 0], [1,0,0], 'r');
            arrow3([0 0 0], [0,1,0], 'b');
            arrow3([0 0 0], [0,0,1], 'g');
        end
    end

end