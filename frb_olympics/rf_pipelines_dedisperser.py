import os
import copy
import numpy as np

import frb_olympics


####################################################################################################


import_successful = False

try:
    import rf_pipelines
    import_successful = True
except ImportError:
    pass


####################################################################################################
#
# frb_olympics_stream:


if import_successful:
    class frb_olympics_stream(rf_pipelines.wi_stream):
        def __init__(self, sparams, nt_chunk=1024):            
            rf_pipelines.wi_stream.__init__(self, 'frb_olympics_stream')
            
            self.sparams = sparams
            self.intensity_arr = None   # see set_intensity() below

            # Note: the following members are properties defined in the C++ base class.
            self.nfreq = sparams.nfreq
            self.nt_chunk = nt_chunk
            

        def set_intensity_arr(self, intensity_arr):
            assert (intensity_arr is None) or intensity_arr.shape == (self.nfreq, self.sparams.nsamples)
            self.intensity_arr = intensity_arr
            

        def _bind_stream(self, json_attrs):
            json_attrs['freq_lo_MHz'] = self.sparams.freq_lo_MHz
            json_attrs['freq_hi_MHz'] = self.sparams.freq_hi_MHz
            json_attrs['dt_sample'] = self.sparams.dt_sample


        def _fill_chunk(self, intensity, weights, pos):
            if self.intensity_arr is None:
                raise RuntimeError('')
                
            assert intensity.shape == (self.nfreq, self.nt_chunk)
            assert weights.shape == (self.nfreq, self.nt_chunk)
            
            # Number of samples to be written
            nt = self.nt_chunk
            nt = min(nt, self.sparams.nsamples - pos)
            nt = max(nt, 0)

            if nt > 0:
                intensity[:,:nt] = self.intensity_arr[:,pos:(pos+nt)]
                weights[:,:nt] = 1.0

            if nt < self.nt_chunk:
                intensity[:,nt:] = 0.0
                weights[:,nt:] = 0.0

            # Reminder: the return value from _fill_chunk() should be True normally,
            # or False if end-of-stream has been reached.
            return (pos + self.nt_chunk < self.sparams.nsamples)

        # Note: jsonize(), from_json() methods are usually necessary to define
        # an rf_pipelines.wi_stream, but we don't need them here (I think!)


####################################################################################################


class rf_pipelines_dedisperser(frb_olympics.dedisperser_base):
    def __init__(self, pipeline, tex_label):
        """The 'pipeline' argument should be an object of class rf_pipelines.pipeline_object."""

        if not import_successful:
            # Rather than throw an exception, we let 'import rf_pipelines' throw an uncaught
            # exception, so that the caller can see what the problem is.

            import rf_pipelines as r
            raise RuntimeError("frb_olympics.rf_pipelines_dedisperser internal error: 'import rf_pipelines' worked on the second try?!")

        if not isinstance(pipeline, rf_pipelines.pipeline_object):
            raise RuntimeError("rf_pipelines_dedisperser.__init__: 'pipeline' constructor argument must be an object of class rf_pipelines.pipeline_object")
        
        if tex_label is None:
            raise RuntimeError("rf_pipelines_dedisperser.__init__: 'tex_label' constructor argument cannot be None")

        # FIXME should add check that pipeline is unbound
        pv = self._analyze_pipeline(pipeline)

        if pv is None:
            raise RuntimeError("rf_pipelines_dedisperser.__init__: pipeline does not appear to contain a dedisperser"
                               + " (this may mean that frb_olympics.rf_pipelines_dedisperser._analyze_pipeline() is out of date)")
        
        frb_olympics.dedisperser_base.__init__(self, tex_label, precomputed_variance=pv)

        self.base_pipeline = pipeline


    def _analyze_pipeline(self, pipeline):
        """
        Helper function called by constructor.  Does some pipeline sanity checking.  
        Return value is:
        
          None: no dedisperser was found in pipeline
          False: dedisperser found, and does not use precomputed analytic variance
          True: dedisperser found, and does use precomputed analytic variance
        """

        if isinstance(pipeline, rf_pipelines.pipeline_object):
            #
            # This is sort of a hack, but we replace the pipeline_object by its jsonization, and use
            # json data structures throughout this routine instead of pipeline_objects.  This is because
            # the container classes in rf_pipelines (e.g. rf_pipelines.pipeline) do not currently define
            # a python API for retreiving their contents, so the only way to "see" inside is by jsonizing.
            #
            # FIXME: when the rf_pipelines python API is more developed, it should be possible to remove
            # this hack.  (One disadvantage of the hack is that all transforms in the pipeline must define
            # jsonize().)

            pipeline = pipeline.jsonize()  # fall through...

            
        if isinstance(pipeline, list):
            ret = None
            count = 0
            
            for p in pipeline:
                t = self._analyze_pipeline(p)
                if t is not None:
                    ret = t
                    count += 1

            if count > 1:
                raise RuntimeError("frb_olympics.rf_pipelines_dedisperser.__init__: pipeline defines multiple dedispersers?!")
            
            return ret

        assert isinstance(pipeline, dict)
        assert pipeline.has_key('class_name')

        if pipeline['class_name'] == 'bonsai_dedisperser_python':
            if not pipeline['track_global_max']:
                raise RuntimeError("rf_pipelines_dedisperser.__init__: 'track_global_max' flag is not set"
                                   + " (this may mean that frb_olympics.rf_pipelines_dedisperser._analyze_pipeline() is out of date)")            
            return pipeline['use_analytic_normalization']

        if pipeline['class_name'] == 'bonsai_dedisperser_cpp':
            raise RuntimeError("rf_pipelines_dedisperser.__init__: pipeline contains a bonsai_dedisperser_cpp, not a bonsai_dedisperser_python"
                               + " (currently, rf_pipelines defines two bonsai dedisperser classes, and only the python class will work in frb_olympics)")

        if pipeline['class_name'] == 'pipeline':
            return self._analyze_pipeline(pipeline['elements'])

        if pipeline['class_name'] == 'wi_sub_pipeline':
            return self._analyze_pipeline(pipeline['sub_pipeline'])

        good_class_names = [ 'badchannel_mask',
                             'intensity_clipper',
                             'mask_expander',
                             'mask_filler',
                             'noise_filler',
                             'pipeline_fork',
                             'polynomial_detrender',
                             'spline_detrender',
                             'std_dev_clipper' ]
        
        if pipeline['class_name'] not in good_class_names:
            print >>sys.stderr, "frb_olympics.rf_pipelines_dedisperser.__init__: unrecognized pipeline_object class '%s' in pipeline" % pipeline['class_name']
            print >>sys.stderr, "   (This may mean that frb_olympics.rf_pipelines_dedisperser._analyze_pipeline() is out of date)"

        return None
        

    def init_search_params(self, sparams):
        """Overrides dedisperser_base.init_search_params().  The 'sparams' argument is an instance of 'class search_params'."""

        if hasattr(self, 'search_params'):
            raise RuntimeError('double call to frb_olympics.bonsai_dedisperser.init_search_params()')

        # FIXME I would like to add sanity checks here on the following search_params:
        #   nfreq, freq_lo_MHz, freq_hi_MHz, dt_sample, max_dm.

        s = frb_olympics_stream(sparams)
        p = rf_pipelines.pipeline([s, self.base_pipeline])
        p.bind(outdir=None, verbosity=0)
        
        self.search_params = sparams
        self.full_pipeline = p
        self.stream = s
        
    
    def allocate(self):
        """Overrides dedisperser_base.allocate()."""
        self.full_pipeline.allocate()


    def dedisperse(self, arr):
        """
        Overrides dedisperser_base.dedisperse().  The 'arr' argument is a float32 array of shape (nfreq, nsamples).

        The return value of dedisperse() is a dictionary with 3 members: 'snr', 'dm', and one of { 'tmid', 'tini', 'tfin' }.
        This allows the dedisperser to use one of three possible definitions of the arrival time of an FRB:
           - "initial" arrival time: arrival time at the highest frequency in the band (i.e. least delayed)
           - "final" arrival time: arrival time at the lowest frequency in the band (i.e. most delayed)
           - "middle" arrival time: average of initial and final times (warning: not the arrival time at the central frequency!)
        """

        self.stream.set_intensity_arr(arr)
        j = self.full_pipeline.run(outdir=None, verbosity=0)        
        self.stream.set_intensity_arr(None)

        expected_keys = [ 'frb_global_max_trigger',
                          'frb_global_max_trigger_dm',
                          'frb_global_max_trigger_tfinal' ]

        for k in expected_keys:
            if not j.has_key(k):
                raise RuntimeError("rf_pipelines_dedisperser: pipeline failed to set json attribute '%s' as expected" % k)
        
        return { 'snr': j['frb_global_max_trigger'],
                 'dm': j['frb_global_max_trigger_dm'],
                 'tfin': j['frb_global_max_trigger_tfinal'] }


    def deallocate(self):
        """Overrides dedisperser_base.deallocate()."""
        self.full_pipeline.deallocate()


    def jsonize(self):
        """
        Overrides dedisperser_base.jsonize().  The return value should be a python dictionary which is valid JSON.
        Note: Caller will add additional members 'module_name', 'class_name', 'tex_label' to the dictionary.
        """
        return { "pipeline": self.base_pipeline.jsonize() }


    @staticmethod
    def from_json(j, filename=None):
        """
        This is the "inverse" of jsonize(): it takes a python dictionary which is valid JSON, 
        and returns an rf_pipelines_dedisperser instance.
        """

        expected_keys = [ 'tex_label', 'pipeline' ]
        r = frb_olympics.json_read_helper(j, filename, 'rf_pipelines_dedipserser.from_json()', expected_keys)
        
        p = rf_pipelines.pipeline_object.from_json(r.json['pipeline'])
        return rf_pipelines_dedisperser(p, r.json['tex_label'])


