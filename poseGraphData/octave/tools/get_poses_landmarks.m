function [poses, landmarks] = get_poses_landmarks(g)

poses = [];
landmarks = [];
values = g.idLookup;
keys = fieldnames(values);
for i = 1:length(keys)
  value = values.(keys{i});
  dim = value.dimension;
  offset = value.offset;
  if (dim == 3)
    poses = [poses; offset];
  elseif (dim == 2)
    landmarks = [landmarks; offset];
  end
end
poses = poses +1 ;
landmarks = landmarks + 1;

end
