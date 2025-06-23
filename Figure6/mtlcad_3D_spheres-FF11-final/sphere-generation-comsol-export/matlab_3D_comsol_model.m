% Model exported on Nov 5 2024, 08:16 by COMSOL 4.4.0.150.
function out = model
    import com.comsol.model.*
    import com.comsol.model.util.*
    model = ModelUtil.create('Model');
    model.modelPath('H:\My Drive\QNN_Group\external-collaboration\2024_Ottavio_comsol\sphere3D-mtlcad\mtlcad_3D_150_spheres-final');

    model.modelNode.create('comp1');
    model.geom.create('geom1', 3);
    model.mesh.create('mesh1', 'geom1');
    model.physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');
    
    model.study.create('std1');
    model.study('std1').feature.create('freq', 'Frequency');
    model.study('std1').feature('freq').activate('ewfd', true);
    
    model.geom('geom1').lengthUnit([native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
    model.geom('geom1').scaleUnitValue(true);
    
    model.param.set('c0', '3e8 [m/s]');
%     model.param.set('lambda', '0.01 [um]');
    model.param.set('lambda', '1 [um]');
    model.param.set('freq0', 'c0/lambda');
    
    model.geom('geom1').feature.create('blk1', 'Block');
%     model.geom('geom1').feature('blk1').setIndex('size', '6+100*lambda', 0);
%     model.geom('geom1').feature('blk1').setIndex('size', '6+100*lambda', 1);
%     model.geom('geom1').feature('blk1').setIndex('size', '6+100*lambda', 2);
    model.geom('geom1').feature('blk1').setIndex('size', '6+lambda', 0);
    model.geom('geom1').feature('blk1').setIndex('size', '6+lambda', 1);
    model.geom('geom1').feature('blk1').setIndex('size', '6+lambda', 2);
    model.geom('geom1').feature('blk1').set('base', 'center');
    model.geom('geom1').feature('blk1').setIndex('pos', '1', 0);
    model.geom('geom1').feature('blk1').setIndex('pos', '1', 1);
    model.geom('geom1').feature('blk1').setIndex('pos', '1', 2);
    model.geom('geom1').feature('blk1').set('layertop', 'on');
    model.geom('geom1').feature('blk1').set('layerback', 'on');
    model.geom('geom1').feature('blk1').set('layerfront', 'on');
    model.geom('geom1').feature('blk1').set('layerright', 'on');
    model.geom('geom1').feature('blk1').set('layerleft', 'on');
%     model.geom('geom1').feature('blk1').setIndex('layer', '50*lambda', 0);
    model.geom('geom1').feature('blk1').setIndex('layer', '0.5*lambda', 0);
    model.geom('geom1').run('blk1');
        
    data = readtable('mtlcad_3D_spheres_FF75_final.csv');
    for k1=1:1
        sphereName = sprintf('sph%d', k1);
        model.geom('geom1').feature.create(sphereName, 'Sphere');
        model.geom('geom1').feature(sphereName).set('r', data{k1,4});
        model.geom('geom1').feature(sphereName).setIndex('pos', data{k1,1}, 0);
        model.geom('geom1').feature(sphereName).setIndex('pos', data{k1,2}, 1);
        model.geom('geom1').feature(sphereName).setIndex('pos', data{k1,3}, 2);
    end
    model.geom('geom1').runPre('fin');
    model.geom('geom1').run;

    model.geom('geom1').feature('sph1').set('createselection', 'on');
    model.geom('geom1').selection.create('csel1', 'CumulativeSelection');
    model.geom('geom1').selection('csel1').name('Sphere');
    model.geom('geom1').feature('sph1').set('contributeto', 'csel1');

    for k1=2:size(data,1)
        sphereName = sprintf('sph%d', k1);
        model.geom('geom1').feature.create(sphereName, 'Sphere');
        model.geom('geom1').feature(sphereName).set('r', data{k1,4});
        model.geom('geom1').feature(sphereName).setIndex('pos', data{k1,1}, 0);
        model.geom('geom1').feature(sphereName).setIndex('pos', data{k1,2}, 1);
        model.geom('geom1').feature(sphereName).setIndex('pos', data{k1,3}, 2);
        model.geom('geom1').feature(sphereName).set('createselection', 'on');
        model.geom('geom1').feature(sphereName).set('contributeto', 'csel1');
    end
    
    model.selection.create('sel1', 'Explicit');
    model.selection('sel1').name('Entire');
    model.selection('sel1').all;

    model.selection.create('box1', 'Box');
    model.selection('box1').name('Box_simulate');
    model.selection('box1').set('xmin', '-2');
    model.selection('box1').set('xmax', '4');
    model.selection('box1').set('ymin', '-2');
    model.selection('box1').set('ymax', '4');
    model.selection('box1').set('zmin', '-2');
    model.selection('box1').set('zmax', '4');
    model.selection('box1').set('condition', 'inside');
    
    model.selection.create('dif1', 'Difference');
    model.selection('dif1').name('PML');
    model.selection('dif1').set('add', {'sel1'});
    model.selection('dif1').set('subtract', {'box1'});
    
    model.coordSystem.create('pml1', 'geom1', 'PML');
    model.coordSystem('pml1').selection.named('dif1');
   
    model.material.create('mat1');
    model.material('mat1').selection.all;
    model.material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
    model.material('mat1').propertyGroup('RefractiveIndex').set('n', {'1'});
    model.material('mat1').name('Air');

    model.material.create('mat2');
    model.material('mat2').selection.named('geom1_csel1_dom');
    model.material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
    model.material('mat2').propertyGroup('RefractiveIndex').set('n', {'1.77'});
    model.material('mat2').name('Alumina');
    
    model.physics('ewfd').feature.create('ecd1', 'ExternalCurrentDensity', 3);
    model.physics('ewfd').feature('ecd1').selection.named('geom1_csel1_dom');
    model.physics('ewfd').feature('ecd1').set('Je', {'0' '0' 'j*(2*pi*c_const/lambda)/Z0_const*1[T]'});
    
    model.mesh('mesh1').feature.create('size1', 'Size');
    model.mesh('mesh1').feature('size').set('hauto', '4');
    model.mesh('mesh1').feature('size1').selection.geom('geom1', 3);
    model.mesh('mesh1').feature('size1').selection.named('geom1_csel1_dom');
    model.mesh('mesh1').feature('size1').set('hauto', '2');
    model.mesh('mesh1').feature.create('cr1', 'CornerRefinement');
    model.mesh('mesh1').feature('cr1').selection('boundary').named('geom1_csel1_bnd');
    model.mesh('mesh1').feature.create('dis1', 'Distribution');
    model.mesh('mesh1').feature('dis1').selection.named('dif1');
    model.mesh('mesh1').feature('dis1').set('numelem', '25');
    model.mesh('mesh1').feature.create('ftet1', 'FreeTet');
    model.mesh('mesh1').run;
    
    model.study('std1').feature('freq').set('plist', 'freq0');
    model.study('std1').feature.create('param', 'Parametric');
    model.study('std1').feature('param').setIndex('pname', 'lambda', 0);
%     model.study('std1').feature('param').setIndex('plistarr', '10^{range(log10(0.01),1/100,log10(1))}', 0);
    model.study('std1').feature('param').setIndex('plistarr', '10^{range(log10(1),1/100,log10(10))}', 0);

    % Save the model
    model.save('H:\My Drive\QNN_Group\external-collaboration\2024_Ottavio_comsol\sphere3D-mtlcad\mtlcad_3D_150_spheres-final\mtlcad_3D_model_FF75_wv2.mph');

out = model;

