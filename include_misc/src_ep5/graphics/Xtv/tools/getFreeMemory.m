function memAvailable = getFreeMemory()
%Try to get available memory on Linux and Windows (in Bytes)
%
% We try different approaches until we get a result that makes sense.
% On Windows, there is memory(), but this is still not available in
% the UNIX builds of Matlab.
% 
% If memory() fails, we search for info in /proc/meminfo.
% If that fails, we use the function free, which should print the free memory.
% But note that
%   free_memory < available_memory_before_swapping
% on Linux at least.
%
% I have no idea what will happen on a Mac.
%
% See also
% memory

% blame: Michael Völker

fallback = 4 * 2^30;    % 2^30 Bytes = 1 GiB
memAvailable = 'meh';   % some invalid default

% valid output must have these attributes
attributes = {'nonempty', 'scalar', 'real', 'finite', 'positive'};

try
    % currently only in windows:
    [~,sys] = memory;
    memAvailable = sys.PhysicalMemory.Available;
catch
    try
        % Linux, version #1
        % https://git.kernel.org/cgit/linux/kernel/git/torvalds/linux.git/commit/?id=34e431b0ae398fc54ea69ff85ec700722c9da773
        % output in kiB:
        [~,memAvailable] = unix('awk ''tolower($1) ~ /^memavailable/{print $(NF-1)}'' /proc/meminfo');
        memAvailable = str2double(memAvailable);
        memAvailable = 1024 * memAvailable;
        validateattributes( memAvailable, {'numeric'}, attributes )
    catch
        try
            % Linux, version #2, less accurate, more portable
            % MemFree + Buffers + Cached + SReclaimable - Shmem
            % output in kiB:
            [~,memAvailable] = unix('awk ''tolower($1) ~ /^memfree|^buffers|^cached|^sreclaimable/{s+=$(NF-1)} tolower($1) ~ /shmem:/{s-=$(NF-1)} END{print s}'' /proc/meminfo');
            memAvailable = str2double(memAvailable);
            memAvailable = 1024 * memAvailable;
            validateattributes( memAvailable, {'numeric'}, attributes )
        catch
            % Linux, version #3, even less accurate
            [~,memAvailable] = unix('LANG=C free -b | awk ''/buffers\/cache/ {print $NF}''');
            memAvailable = str2double(memAvailable);
            % note: memAvailable is intentionally unchecked here (could be empty, NaN, ...)
        end
    end
end

try
    validateattributes( memAvailable, {'numeric'}, attributes )
catch
    warning( [mfilename() ':NoLuck'], 'Could not find available memory. Setting to %g MiB...', fallback/2^20 )
    memAvailable = fallback;
end

end % of getFreeMemory()