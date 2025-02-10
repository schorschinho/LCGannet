% coreg_siemens.m
% Georg Oeltzschner, Johns Hopkins University 2019.
% 
% USAGE:
% [vol_mask] = coreg_siemens(in, vol_image, maskFile);
% 
% DESCRIPTION:
% Creates a SPM volume containing a voxel mask with the same dimensions 
% as the SPM volume containing a structural image. The voxel mask will be
% created with the geometry information stored in the FID-A structure "in".
% This routine works for the Siemens TWIX, DICOM, and RDA data types, where
% the voxel orientation is defined by a normal vector and a rotation around
% this normal vector.
% 
% CREDITS:
% The routine for correct determination of the phase and readout
% directions of the MRS voxel is adapted from 
% vox2ras_rsolveAA.m
% (Dr. Rudolph Pienaar, Massachusetts General Hospital, Boston)
%
% INPUTS:
% in        = Input data structure.
% vol_image	= SPM volume of the structural image that the voxel defined in
%           the 'geometry' field of the input data structure 'in' is
%           supposed to be co-registered to.
% maskFile  = Filename under which the SPM volume of the co-registered
%           voxel mask is supposed to be saved.
%
% OUTPUTS:
% vol_mask  = SPM volume of the coregistered voxel mask.
% T1_max    = maximum intensity of the image volume.

function [vol_mask, T1_max, voxel_ctr] = coreg_siemens(in, vol_image, maskFile)

% Deactivate MATLAB warnings and load geometry parameters
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:qhullmx:InternalWarning');
geom = in.geometry;

%%% 1. PREPARE THE STRUCTURAL IMAGE
% Create SPM volume and read in the NIfTI file with the structural image.
[T1,XYZ]    = spm_read_vols(vol_image);
T1_max      = max(T1(:));
%Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
%tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
[~,voxdim] = spm_get_bbox(vol_image,'fv');
voxdim = abs(voxdim)';
halfpixshift = -voxdim(1:3)/2;
halfpixshift(3) = -halfpixshift(3);
XYZ = XYZ + repmat(halfpixshift, [1 size(XYZ,2)]);


%%% 2. GENERATE THE COORDINATES OF THE VOXEL CORNERS
% Extract voxel position and rotation parameters
NormSag         = geom.rot.NormSag;
NormCor         = geom.rot.NormCor;
NormTra         = geom.rot.NormTra;
VoI_InPlaneRot  = geom.rot.VoI_InPlaneRot;
% Correct voxel offsets by table position (if field exists)
if isfield(geom.pos,'TablePosTra')
    VoxOffs = [geom.pos.PosSag geom.pos.PosCor geom.pos.PosTra+geom.pos.TablePosTra];
else
    VoxOffs = [geom.pos.PosSag geom.pos.PosCor geom.pos.PosTra];
end

% Parse direction cosines of the MRS voxel's normal vector and the rotation angle
% around the normal vector
% The direction cosine is the cosine of the angle between the normal
% vector and the respective direction.
% Example: If the normal vector points exactly along the FH direction, then: 
% NormSag = cos(90) = 0, NormCor = cos(90) = 0, NormTra = cos(0) = 1.
Norm = [-NormSag -NormCor NormTra];
ROT = VoI_InPlaneRot;
% Find largest element of normal vector of the voxel to determine primary
% orientation. 
% Example: if NormTra has the smallest out of the three Norm
% values, the angle of the normal vector with the Tra direction (FH) is the
% smallest, and the primary orientation is transversal.
[~, maxdir] = max([abs(NormSag) abs(NormCor) abs(NormTra)]);
switch maxdir
    case 1
        vox_orient = 's'; % 't' = transversal, 's' = sagittal', 'c' = coronal;
    case 2
        vox_orient = 'c'; % 't' = transversal, 's' = sagittal', 'c' = coronal;
    case 3
        vox_orient = 't'; % 't' = transversal, 's' = sagittal', 'c' = coronal;
end
    
% Phase reference vector
% Adapted from Rudolph Pienaar's "vox2ras_rsolveAA.m" and
% Andre van der Kouwe's "autoaligncorrect.cpp"
Phase	= zeros(3, 1);
switch vox_orient
    case 't'
        % For transversal voxel orientation, the phase reference vector lies in
        % the sagittal plane
        Phase(1)	= 0;
        Phase(2)	=  Norm(3)*sqrt(1/(Norm(2)*Norm(2)+Norm(3)*Norm(3)));
        Phase(3)	= -Norm(2)*sqrt(1/(Norm(2)*Norm(2)+Norm(3)*Norm(3)));
        VoxDims = [geom.size.VoI_PeFOV geom.size.VoI_RoFOV geom.size.VoIThickness];
    case 'c'
        % For coronal voxel orientation, the phase reference vector lies in
        % the transversal plane
        Phase(1)	=  Norm(2)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(2)	= -Norm(1)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(3)	= 0;
        VoxDims = [geom.size.VoI_PeFOV geom.size.VoI_RoFOV geom.size.VoIThickness];
    case 's'
        % For sagittal voxel orientation, the phase reference vector lies in
        % the transversal plane
        Phase(1)	= -Norm(2)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(2)	=  Norm(1)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(3)	= 0;
        VoxDims = [geom.size.VoI_PeFOV geom.size.VoI_RoFOV geom.size.VoIThickness];
end

% The readout reference vector is the cross product of Norm and Phase
Readout = cross(Norm, Phase);
M_R = zeros(4, 4);
M_R(1:3, 1)	= Phase;
M_R(1:3, 2)	= Readout;
M_R(1:3, 3) = Norm;

% Define matrix for rotation around in-plane rotation angle
M3_Mu	= [	 cos(ROT)	sin(ROT)	0
            -sin(ROT)	cos(ROT)	0
            0           0           1];
        
M3_R	= M_R(1:3,1:3)	* M3_Mu;
M_R(1:3,1:3)	= M3_R;

% The MGH vox2ras matrix inverts the Readout column
M_R		= M_R *   [ 1  0  0  0
                    0 -1  0  0
                    0  0  1  0
                    0  0  0  1];

% Final rotation matrix
rotmat = M_R(1:3,1:3);

% We need to flip ap and lr axes to match NIFTI convention
VoxOffs(1) = -VoxOffs(1);
VoxOffs(2) = -VoxOffs(2);

% Define voxel coordinates before rotation and transition
vox_ctr = ...
    [VoxDims(1)/2 -VoxDims(2)/2  VoxDims(3)/2;
    -VoxDims(1)/2 -VoxDims(2)/2  VoxDims(3)/2;
    -VoxDims(1)/2  VoxDims(2)/2  VoxDims(3)/2;
     VoxDims(1)/2  VoxDims(2)/2  VoxDims(3)/2;
    -VoxDims(1)/2  VoxDims(2)/2 -VoxDims(3)/2;
     VoxDims(1)/2  VoxDims(2)/2 -VoxDims(3)/2;
     VoxDims(1)/2 -VoxDims(2)/2 -VoxDims(3)/2;
    -VoxDims(1)/2 -VoxDims(2)/2 -VoxDims(3)/2];

% Apply rotation as prescribed
vox_rot = rotmat*vox_ctr.';

% Shift rotated voxel by the center offset to its final position
vox_ctr_coor = [VoxOffs(1) VoxOffs(2) VoxOffs(3)];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
vox_corner = vox_rot + vox_ctr_coor;


%%% 3. CREATE AND SAVE THE VOXEL MASK
% Create a mask with all voxels that are inside the voxel
mask = zeros(1,size(XYZ,2));
sphere_radius = sqrt((VoxDims(1)/2)^2+(VoxDims(2)/2)^2+(VoxDims(3)/2)^2);
distance2voxctr = sqrt(sum((XYZ-repmat([VoxOffs(1) VoxOffs(2) VoxOffs(3)].',[1 size(XYZ, 2)])).^2,1));
sphere_mask(distance2voxctr <= sphere_radius) = 1;
mask(sphere_mask == 1) = 1;
XYZ_sphere = XYZ(:,sphere_mask == 1);
tri = delaunayn([vox_corner.'; [VoxOffs(1) VoxOffs(2) VoxOffs(3)]]);
tn = tsearchn([vox_corner.'; [VoxOffs(1) VoxOffs(2) VoxOffs(3)]], tri, XYZ_sphere.');
isinside = ~isnan(tn);
mask(sphere_mask==1) = isinside;

% Take over the voxel dimensions from the structural
mask = reshape(mask, vol_image.dim);

% Fill in the SPM volume header information
vol_mask.fname   = maskFile;
vol_mask.dim     = vol_image.dim;
vol_mask.dt      = vol_image.dt;
vol_mask.mat     = vol_image.mat;
vol_mask.descrip = 'MRS_voxel_mask';

% Write the SPM volume to disk
vol_mask = spm_write_vol(vol_mask,mask);

% Store voxel centre for output figure
voxel_ctr = [VoxOffs(1) VoxOffs(2) VoxOffs(3)];

% Reactivate MATLAB warnings
warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:qhullmx:InternalWarning');

end
