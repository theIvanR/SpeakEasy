function geom = load_mesh(mesh_input, scale)
%LOAD_MESH  Convert STL / triangulation / face-vertex data into geom.
%
% Supports:
%   - triangulation object
%   - struct with fields faces, vertices
%   - .mat file containing geom, or faces+vertices
%   - .stl file if stlread is available

    if nargin < 2 || isempty(scale)
        scale = 1;
    end

    if isa(mesh_input, 'triangulation')
        faces = mesh_input.ConnectivityList;
        vertices = mesh_input.Points * scale;
        geom = mesh_to_geom(vertices, faces);
        return;
    end

    if isstruct(mesh_input) && isfield(mesh_input, 'faces') && isfield(mesh_input, 'vertices')
        faces = mesh_input.faces;
        vertices = mesh_input.vertices * scale;
        geom = mesh_to_geom(vertices, faces);
        return;
    end

    if isstring(mesh_input) || ischar(mesh_input)
        [~, ~, ext] = fileparts(char(mesh_input));

        switch lower(ext)
            case '.mat'
                S = load(mesh_input);
                if isfield(S, 'geom')
                    geom = S.geom;
                    return;
                elseif isfield(S, 'faces') && isfield(S, 'vertices')
                    geom = mesh_to_geom(S.vertices * scale, S.faces);
                    return;
                else
                    error('MAT file must contain geom or faces+vertices');
                end

            case '.stl'
                assert(exist('stlread', 'file') == 2, 'stlread not found on path');
                TR = stlread(mesh_input);
                if isa(TR, 'triangulation')
                    faces = TR.ConnectivityList;
                    vertices = TR.Points * scale;
                    geom = mesh_to_geom(vertices, faces);
                elseif isstruct(TR) && isfield(TR, 'ConnectivityList') && isfield(TR, 'Points')
                    geom = mesh_to_geom(TR.Points * scale, TR.ConnectivityList);
                else
                    error('Unsupported stlread output');
                end
                return;

            otherwise
                error('Unsupported mesh file type: %s', ext);
        end
    end

    error('Unsupported mesh_input type');
end

function geom = mesh_to_geom(vertices, faces)
%MESH_TO_GEOM  Build centroid/area geometry from triangular faces.

    nF = size(faces, 1);
    xc = zeros(nF, 3);
    A  = zeros(nF, 1);
    aeq = zeros(nF, 1);

    for i = 1:nF
        tri = vertices(faces(i, :), :);
        xc(i, :) = mean(tri, 1);

        v1 = tri(2, :) - tri(1, :);
        v2 = tri(3, :) - tri(1, :);
        Ai = 0.5 * norm(cross(v1, v2));

        A(i) = Ai;
        aeq(i) = sqrt(Ai / pi);
    end

    geom = struct();
    geom.xc  = xc;
    geom.A   = A;
    geom.aeq = aeq;
    geom.Np  = nF;
    geom.faces = faces;
    geom.vertices = vertices;
end