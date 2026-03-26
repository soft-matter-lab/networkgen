function Nvec = AssignPerBondMono(obj, Bonds, Atoms)
% -------------------------------------------------------------------------
% AssignPerBondMono
% - Assign one identical Kuhn segment count to all bonds
%
% Uses:
%   obj.perbond.kuhn.mono.value
%
% OUTPUT:
%   Nvec : [Nbonds x 1]
% -------------------------------------------------------------------------

    nbonds = size(Bonds,1);
    natoms = size(Atoms,1);

    Nvec = zeros(nbonds,1);
    if nbonds == 0
        return;
    end

    Nmono = obj.perbond.kuhn.mono.value;

    if isempty(Nmono)
        Nmono = 20;
    end

    Nmono = round(Nmono);
    Nmono = max(Nmono, 1);

    Nvec(:) = Nmono;

    obj.log.print('   Kuhn assignment mode: mono\n');
    obj.log.print('   Kuhn-to-crosslinker ratio %0.4f\n', sum(Nvec)/max(natoms,1));
    obj.log.print('   Average chain length %0.4f\n', sum(Nvec)/max(length(Nvec),1));

end