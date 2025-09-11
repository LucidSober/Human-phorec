function status = record_dt(t,~,flag)
  persistent lastT
  if strcmp(flag,'init'), lastT = []; setappdata(0,'dt_prev',NaN);
  elseif isempty(flag)    % 每个"已接受"步才会来这里
    if ~isempty(lastT)
      setappdata(0,'dt_prev', t(end)-lastT);
    end
    lastT = t(end);
  end
  status = 0;
end