% Identity test for handle objects (Octave crashes on isequal() for cyclic handle graphs).
% Objects get a lazily assigned unique id stored in their 'uid' property.
function tf=isSame(a,b)
    tf=false;
    if isempty(a)||isempty(b)||~isobject(a)||~isobject(b)
        return
    end
    tf=(objUid(a)==objUid(b));
end

function id=objUid(h)
    persistent counter
    if isempty(counter)
        counter=0;
    end
    if isempty(h.uid)
        counter=counter+1;
        h.uid=counter;
    end
    id=h.uid;
end
