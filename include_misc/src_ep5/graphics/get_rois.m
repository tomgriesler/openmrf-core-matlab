function [rois ]= get_rois(varargin)

parseVarargin

  [h]=gui_test_2(varargin{:});
  
  f=guidata(h);
  rois = f.rois;
