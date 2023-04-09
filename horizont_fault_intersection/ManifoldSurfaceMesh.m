classdef ManifoldSurfaceMesh < handle
    properties (Constant)
        INVALID_IND = intmax('uint64'); % max of size_t
    end
    properties
        m_nfaces
        m_nvertices
        m_nhed_count = uint64(0);
        m_nedge_count = uint64(0);
        m_nboundary_loops = uint64(0);

        m_vHalfEdgeArr % returns the halfedge id for a vertex
        m_fHalfEdgeArr % returns the halfedge id for a face

        %Resizable arrays in the original code, here these array are
        %automatically resizable
        m_heNextArr = [uint64(1)]; % return the hed next id for a halfedge
        m_heVertexArr = [uint64(1)]; % return the first vertice of a hed
        m_heFaceArr = [uint64(1)]; % return the face from where the hed belogs
        % interior heds
        m_ninteriorHeds
        %exterior_heds = m_nhed_count - m_ninteriorHeds
    end
    methods
        function this = ManifoldSurfaceMesh(faces)
           
            this.m_nfaces = uint64(size(faces,1));
            this.m_nvertices = uint64(max(max(faces)));

            this.m_vHalfEdgeArr = ManifoldSurfaceMesh.INVALID_IND * uint64(ones(this.m_nvertices, 1)); % get a handle for any halfedge
            this.m_fHalfEdgeArr = ManifoldSurfaceMesh.INVALID_IND * uint64(ones(this.m_nfaces, 1)); %get a handle for any halfedge

            hed_hash = containers.Map('KeyType','char','ValueType','uint64');

            for id_face = 1:this.m_nfaces
                face_degree = 3; %just triangular faces
                first_hed_id = ManifoldSurfaceMesh.INVALID_IND;
                prev_hed_id = ManifoldSurfaceMesh.INVALID_IND;

                for id_face_hed = 1:face_degree
                    id_v0 = faces(id_face, id_face_hed);
                    id_v1 = faces(id_face, mod(id_face_hed,face_degree) + 1);

                    hedKey = [id_v0, id_v1];
                    hedTwinKey = [id_v1, id_v0];

                    hed_id = this.createHedId(hed_hash, hedKey);
                    hed_id_twin = this.createHedId(hed_hash, hedTwinKey);

                    %check if twin is already created. Only the twin hed is
                    %checked because, the other hed we already know that
                    %exists. 
                    if (hed_id_twin == ManifoldSurfaceMesh.INVALID_IND)
                        %Only the original hed is creating a hash, the
                        %twin hash is created only when we find twice
                        hed_id = this.updateHedCount();
                        hed_hash(this.key_enc(hedKey)) = hed_id; % this would be simple if hed_id was a rvalue

                        hed_id_twin = this.getTwin(hed_id);
                        this.m_heNextArr(hed_id) = ManifoldSurfaceMesh.INVALID_IND;
                        this.m_heNextArr(hed_id_twin) = ManifoldSurfaceMesh.INVALID_IND;
                        
                        this.m_heVertexArr(hed_id) = id_v0;
                        this.m_heVertexArr(hed_id_twin) = id_v1;
                        
                        %We init the face array for the twin also, but at
                        %the boundary loop we check if it is still invalid
                        this.m_heFaceArr(hed_id) = ManifoldSurfaceMesh.INVALID_IND;
                        this.m_heFaceArr(hed_id_twin) = ManifoldSurfaceMesh.INVALID_IND;
                    else
                        %Only in this moment is the twin hed giving a key
                        hed_id= this.getTwin(hed_id_twin);
                        hed_hash(this.key_enc(hedKey)) = hed_id;
                    end

                    this.m_heFaceArr(hed_id) = id_face;

                    % we just need one pointer for the vertex
                    this.m_vHalfEdgeArr(id_v0) = hed_id;

                    if (id_face_hed == 1)
                        %we just need one pointer for the face
                        this.m_fHalfEdgeArr(id_face) = hed_id; 
                        %save first hed to update its next
                        first_hed_id = hed_id;
                    else
                        this.m_heNextArr(prev_hed_id) = hed_id;                        
                    end                    
                    prev_hed_id = hed_id;
                end
                %It remains to point the next hed of the last face hed (its
                %next hed is the first one, due to this we save the
                %first_hed_id
                this.m_heNextArr(prev_hed_id) = first_hed_id;
            end

            %% Treat Boundary Loops
            this.m_ninteriorHeds = this.m_nhed_count;
            for hed_id = uint64(1:this.m_nhed_count)                
                if (this.m_heFaceArr(hed_id) ~= ManifoldSurfaceMesh.INVALID_IND)
                    % heds that belongs to a face
                    continue;
                end
                %now we have a boundary loop
                %note that the boundary loop is also face that comes after
                %all other faces
                this.m_nboundary_loops = this.m_nboundary_loops + 1;                
                boundary_loop_id = this.m_nfaces + this.m_nboundary_loops;
                %consider this hed as the first hed of this "boundary face"
                this.m_fHalfEdgeArr(boundary_loop_id) = hed_id;

                %init walk around to fix the m_heFaceArr (all Invalid faces
                %will change to boundary_loop_id)
                curr_hed = hed_id;
                prev_hed = ManifoldSurfaceMesh.INVALID_IND;
                count_heds = 0;

                while true
                     count_heds = count_heds + 1;

                    %change face boundary hed to loop face
                    this.m_heFaceArr(curr_hed) = boundary_loop_id;

                    %A interesting condition is that all vertices has its
                    %first hed as a interior hed, so we will enforce this
                    %with the twin
                    curr_hed_twin = this.getTwin(curr_hed);
                    vertex_id = this.m_heVertexArr(curr_hed_twin);
                    %enforce the twin
                    this.m_vHalfEdgeArr(vertex_id) = curr_hed_twin;
                    
                    %update interior heds
                    this.m_ninteriorHeds  = this.m_ninteriorHeds  - 1;
                    
                    %Init next operations for boundary loops
                    prev_hed = curr_hed;
                    curr_hed = this.getTwin(this.m_heNextArr(curr_hed_twin));
                    inter_loop = 0;
                    while (this.m_heFaceArr(curr_hed) ~= ManifoldSurfaceMesh.INVALID_IND)                        
                        curr_hed = this.getTwin(this.m_heNextArr(curr_hed));
                        %check if is already the end of the loop
                        if (curr_hed == hed_id)
                            break;
                        end
                        inter_loop = inter_loop + 1;
                        assert (inter_loop < this.m_nhed_count, 'boundary infinite loop');
                    end
                    this.m_heNextArr(curr_hed) = prev_hed;
                    %check if is already the end of the loop
                    if (curr_hed == hed_id)
                        break;
                    end
                    assert (count_heds < this.m_nhed_count, 'boundary infinite loop');
                end
                count_heds;
            end
        end
        function res = updateHedCount(this)
            this.m_nhed_count = this.m_nhed_count + 2;
            this.m_nedge_count = this.m_nedge_count + 1;
            res = uint64(this.m_nhed_count - 1);
        end
        function adjacent_faces = get_adjacentFaces(this,face_id,plotter)
            first_hed = this.m_fHalfEdgeArr(face_id);
            v0 = this.m_heVertexArr(first_hed);
            next_hed = this.m_heNextArr(first_hed);
            v1 = this.m_heVertexArr(next_hed);
            next_hed = this.m_heNextArr(next_hed);
            v2 = this.m_heVertexArr(next_hed);
            star_0 = this.getStar(v0);
            star_1 = this.getStar(v1);
            %plotter.showMesh('blue',star_1);

            star_2 = this.getStar(v2);
            adjacent_faces = unique([star_0, star_1, star_2]);
        end
        function star = getStar(this, node_id)
            first_hed = this.m_vHalfEdgeArr(node_id);            
            star = this.m_heFaceArr(first_hed);
            curr_hed = this.getTwin(this.m_heNextArr(this.m_heNextArr(first_hed)));
            while (curr_hed ~= first_hed)
                face_id = this.m_heFaceArr(curr_hed);
                if face_id > this.m_nfaces %boundary vertex
                    break;
                end
                star = [star, this.m_heFaceArr(curr_hed)];
                curr_hed = this.getTwin(this.m_heNextArr(this.m_heNextArr(curr_hed)));
            end
        end
        function star = getStarVertices(this, node_id)
            first_hed = this.m_vHalfEdgeArr(node_id);            
            next_hed = this.m_heNextArr(this.m_heNextArr(first_hed));
            star = [];
            curr_hed = this.getTwin(next_hed);
            is_boundary = false;
            while (curr_hed ~= first_hed)
                face_id = this.m_heFaceArr(curr_hed);
                star = [star, this.m_heVertexArr(next_hed)];
                next_hed = this.m_heNextArr(this.m_heNextArr(curr_hed));
                curr_hed = this.getTwin(next_hed);
                if face_id > this.m_nfaces %boundary vertex
                    is_boundary = true;
                    break;
                end
            end
            % check the other half
            if (is_boundary)
                new_hed = this.getTwin(first_hed);
                star = [star, this.m_heVertexArr(new_hed)];
                face_id = this.m_heFaceArr(new_hed);
                next_hed = (this.m_heNextArr(new_hed));
                curr_hed = this.getTwin(next_hed);
                while (face_id < this.m_nfaces)
                    face_id = this.m_heFaceArr(curr_hed);
                    star = [star, this.m_heVertexArr(curr_hed)];
                    next_hed = this.m_heNextArr(curr_hed);
                    curr_hed = this.getTwin(next_hed);
                end
            end
        end
        function [bool, v_id] = isJoinedByVertex(this,tri_id_1, tri_id_2)
            bool = false;
            face_degree = 3;
            vetices_id = uint64(zeros(1,face_degree));
            hed = this.m_fHalfEdgeArr(tri_id_1);
            for i = 1:face_degree
                v_id = this.m_heVertexArr(hed);
                vertices_id(i) = v_id;
                hed = this.m_heNextArr(hed);
            end
            hed = this.m_fHalfEdgeArr(tri_id_2);
            for i = 1:face_degree
                v_id = this.m_heVertexArr(hed);
                if ismember(v_id , vertices_id)
                    bool = true;
                    break;
                end
                hed = this.m_heNextArr(hed);
            end
        end
        function [bool, inc] = isJoinedByEdge(this,tri_id_1, tri_id_2)
            bool = false;
            inc = [-1,-1];
            face_degree = 3;
            hed = this.m_fHalfEdgeArr(tri_id_1);            
            for i = 1:face_degree
                twin = this.getTwin(hed);
                face_id = this.m_heFaceArr(twin);
                if (face_id == tri_id_2)
                    bool = true;
                    inc = [this.m_heVertexArr(hed), this.m_heVertexArr(twin)];
                    break;
                end
                hed = this.m_heNextArr(hed);
            end
        end
        function faces_id = getFacesFromBoundary(this, boundary_id)
            %Run BFS
            first_hed = this.m_fHalfEdgeArr(boundary_id);
            face_mark = ~logical(1:this.m_nfaces);
            vertices_mark = ~logical(1:this.m_nvertices);
            first_node = this.m_heVertexArr(first_hed);
            Q = first_node;
            while (~isempty(Q))
                u = Q(1);
                vertices_mark(u) = true;
                Q(1) = [];
                children = this.getStarVertices(u);
                adjacent_faces = this.getStar(u);
                face_mark(adjacent_faces) = true;
                for child = children
                    if (~vertices_mark(child))
                        Q(end+1) = child;
                    end
                end
            end
            range = 1:this.m_nfaces;
            faces_id = range(face_mark);            
        end
    end
    methods (Static)
        function res = createHedId(hash, key)
            key = ManifoldSurfaceMesh.key_enc(key);
            if ~isKey(hash, key)
                hash(key) = ManifoldSurfaceMesh.INVALID_IND;
            end
            res = (hash(key));
        end
        % Function to transform halfedges incidence in chars (works only
        % with row vectors
        % see https://www.mathworks.com/matlabcentral/answers/553150-dictionary-map-with-array-as-key
        function key1_char =  key_enc(key1)
            [~,n] = size(key1);
            assert(n == 2);
            key1_char = ['0x',dec2hex(key1(1)),'0x',dec2hex(key1(2))];
        end
        function key1_dbl =  key_dec(key1_char)
            n = numel(key1_char);
            key1_dbl = transpose(hex2num(transpose(reshape(key1_char,16,n/16))));
        end
        %Bitwise operation to get hed twin
        function hed_twin = getTwin(hed_id)
            hed_twin = bitxor(hed_id - 1, 1) + 1;
        end
    end
end
