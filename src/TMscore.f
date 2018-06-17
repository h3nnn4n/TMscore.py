*************************************************************************
*     This program is to compare two protein structures and identify the
*     best superposition that has the highest TM-score. Input structures
*     must be in the PDB format. By default, TM-score is normalized by
*     the second protein. Users can obtain a brief instruction by simply
*     running the program without arguments. For comments/suggestions,
*     please contact email: zhng@umich.edu.
*
*     Reference:
*     Yang Zhang, Jeffrey Skolnick, Proteins, 2004 57:702-10.
*
*     Permission to use, copy, modify, and distribute this program for
*     any purpose, with or without fee, is hereby granted, provided that
*     the notices on the head, the reference information, and this
*     copyright notice appear in all copies or substantial portions of
*     the Software. It is provided "as is" without express or implied
*     warranty.
******************* Updating history ************************************
*     2005/10/19: the program was reformed so that the score values.
*                 are not dependent on the specific compilers.
*     2006/06/20: selected 'A' if there is altLoc when reading PDB file.
*     2007/02/05: fixed a bug with length<15 in TMscore_32.
*     2007/02/27: rotation matrix from Chain-1 to Chain-2 was added.
*     2007/12/06: GDT-HA score was added, fixed a bug for reading PDB.
*     2010/08/02: A new RMSD matrix was used and obsolete statement removed.
*     2011/01/03: The length of pdb file names were extended to 500.
*     2011/01/30: An open source license is attached to the program.
*     2012/05/07: Improved RMSD calculation subroutine which speeds up
*                 TM-score program by 30%.
*     2012/06/05: Added option '-l L' which calculates TM-score (and maxsub
*                 and GDT scores) normalized by a specific length 'L'.
*     2012/12/17: Added 'TM.sup_atm' to superpose full-atom structures.
*                 The former superposition is for CA-trace only.
*     2013/05/08: Update TM-score so that it can read all alternate location
*                 indicators and residue insertions.
*     2013/05/11: Fix a bug in array overflow.
*     2016/03/23: Extended the program to allow calculating TM-score for for
*                 complex structure comparisons, where multiple-chains are
*                 merged into a single chain. Chain ID is now included in
*                 the output files.
*************************************************************************

c        1         2         3         4         5         6         7 !
c 3456789012345678901234567890123456789012345678901234567890123456789012345678

      program TMscore
      PARAMETER(nmax=5000)

      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0,d0_fix
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)

      character*500 fnam,pdb(100),outname
      character*3 aa(-1:20),seqA(nmax),seqB(nmax)
      character*500 s,du
      character*1 chA(nmax),chB(nmax)
      character*1 chain1(90000),chain2(90000)
      character seq1A(nmax),seq1B(nmax),ali(nmax),du1
      character sequenceA(nmax),sequenceB(nmax),sequenceM(nmax)
      character ins1(nmax),ins2(nmax),ains1(90000),ains2(90000)

      dimension L_ini(100),iq(nmax)
      common/scores/score,score_maxsub,score_fix,score10
      common/GDT/n_GDT05,n_GDT1,n_GDT2,n_GDT4,n_GDT8
      double precision score,score_max,score_fix,score_fix_max
      double precision score_maxsub,score10
      dimension xa(nmax),ya(nmax),za(nmax)

      character*10 aa1,ra1,aa2,ra2
      dimension ia1(90000),aa1(90000),ra1(90000),ir1(90000)
      dimension xa1(90000),ya1(90000),za1(90000)
      dimension ia2(90000),aa2(90000),ra2(90000),ir2(90000)
      dimension xa2(90000),ya2(90000),za2(90000)

      dimension nres1(nmax,32:122),nres2(nmax,32:122) !number of atoms
      character*500 atom1(nmax,30),atom2(nmax,30) !atom name
c     here atom1(i,j) should be atom2(i,j,k) but will beyond dimension

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc

      data aa/ 'BCK','GLY','ALA','SER','CYS',
     &     'VAL','THR','ILE','PRO','MET',
     &     'ASP','ASN','LEU','LYS','GLU',
     &     'GLN','ARG','HIS','PHE','TYR',
     &     'TRP','CYX'/
      character*1 slc(-1:20)
      data slc/'X','G','A','S','C',
     &     'V','T','I','P','M',
     &     'D','N','L','K','E',
     &     'Q','R','H','F','Y',
     &     'W','C'/

*****instructions ----------------->
      call getarg(1,fnam)
      if(fnam.eq.' '.or.fnam.eq.'?'.or.fnam.eq.'-h')then
         write(*,*)
         write(*,*)'Brief instruction for running TM-score program:'
         write(*,*)'(For detail: Zhang & Skolnick,  Proteins, 2004',
     &        ' 57:702-10)'
         write(*,*)
         write(*,*)'1. Run TM-score to compare ''model'' and ',
     &        '''native'':'
         write(*,*)'   >TMscore model native'
         write(*,*)
         write(*,*)'2. Run TM-score to compare two complex structures',
     &        ' with multiple chains'
         write(*,*)'   >TMscore -c model native'
         write(*,*)
         write(*,*)'2. TM-score normalized with an assigned scale d0',
     &        ' e.g. 5 A:'
         write(*,*)'   >TMscore model native -d 5'
         write(*,*)
         write(*,*)'3. TM-score normalized by a specific length, ',
     &        'e.g. 120 AA:'
         write(*,*)'   >TMscore model native -l 120'
         write(*,*)
         write(*,*)'4. TM-score with superposition output, e.g. ',
     &        '''TM.sup'':'
         write(*,*)'   >TMscore model native -o TM.sup'
         write(*,*)'   To view the superimposed CA-traces by rasmol:'
         write(*,*)'      >rasmol -script TM.sup'
         write(*,*)'   To view superimposed atomic models by rasmol:'
         write(*,*)'      >rasmol -script TM.sup_atm'
         write(*,*)
         goto 9999
      endif

******* options ----------->
      m_out=-1
      m_fix=-1
      m_len=-1
      narg=iargc()
      i=0
      j=0
      m_complex=0
 115  continue
      i=i+1
      call getarg(i,fnam)
      if(fnam.eq.'-o')then
         m_out=1
         i=i+1
         call getarg(i,outname)
      elseif(fnam.eq.'-d')then
         m_fix=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)d0_fix
      elseif(fnam.eq.'-c')then
         m_complex=1
      elseif(fnam.eq.'-l')then
         m_len=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)l0_fix
      else
         j=j+1
         pdb(j)=fnam
      endif
      if(i.lt.narg)goto 115

ccccccccc read data from first CA file:

      open(unit=10,file=pdb(1),status='old')
      i=0
      na1=0
      ntt=0
 101  read(10,104,end=102) s
      if(m_complex.eq.0)then
         if(s(1:3).eq.'TER') goto 102
      endif
      if(s(1:4).eq.'ATOM')then
***
         if(s(27:27).eq.' ')then
            du1=' '
         else
            read(s(27:27),*)du1 !Residue_Insertion_tag
         endif
         mk=1
         if(s(17:17).ne.' ')then !ignor alternative residues/atoms
            read(s(13:16),*)du
            read(s(23:26),*)i8
            do i1=1,nres1(i8,ichar(du1))
               if(du.eq.atom1(i8,i1))then !atom is a ALTloc
                  mk=-1
               endif
            enddo
         endif
         if(mk.eq.1)then
            read(s(13:16),*)du
            read(s(23:26),*)i8
            if(nres1(i8,ichar(du1)).lt.30)then
               nres1(i8,ichar(du1))=nres1(i8,ichar(du1))+1
               atom1(i8,nres1(i8,ichar(du1)))=du
            endif
***
            ntt=ntt+1
            if(ntt.ge.90000)goto 102
            na1=na1+1
            read(s,8999)du,ia1(na1),du,aa1(na1),du,ra1(na1),du,
     &           chain1(na1),ir1(na1),du,xa1(na1),ya1(na1),za1(na1)
            ains1(na1)=du1
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.s(13:16).
     &           eq.'  CA')then
               i=i+1
               read(s,103)du,seqA(i),du,chA(i),nresA(i),du,
     &              xa(i),ya(i),za(i)
               ins1(i)=du1
               do j=-1,20
                  if(seqA(i).eq.aa(j))then
                     seq1A(i)=slc(j)
                     goto 21
                  endif
               enddo
               seq1A(i)=slc(-1)
 21            continue
               if(i.ge.nmax)goto 102
            endif
         endif
      endif
      goto 101
 102  continue
 103  format(A17,A3,A1,A1,i4,A4,3F8.3)
 104  format(A100)
 8999 format(a6,I5,a1,A4,a1,A3,a1,a1,I4,a4,3F8.3)
      close(10)
      nseqA=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ccccccccc read data from first CA file:
      open(unit=10,file=pdb(2),status='old')
      i=0
      na2=0
      ntt=0
 201  read(10,204,end=202) s
      if(m_complex.eq.0)then
         if(s(1:3).eq.'TER') goto 202
      endif
      if(s(1:4).eq.'ATOM')then
***
         if(s(27:27).eq.' ')then
            du1=' '
         else
            read(s(27:27),*)du1
         endif
         mk=1
         if(s(17:17).ne.' ')then
            read(s(13:16),*)du
            read(s(23:26),*)i8
            do i1=1,nres2(i8,ichar(du1))
               if(du.eq.atom2(i8,i1))then
                  mk=-1
               endif
            enddo
         endif
         if(mk.eq.1)then
            read(s(13:16),*)du
            read(s(23:26),*)i8
            if(nres2(i8,ichar(du1)).lt.30)then
               nres2(i8,ichar(du1))=nres2(i8,ichar(du1))+1
               atom2(i8,nres2(i8,ichar(du1)))=du
            endif
***
            ntt=ntt+1
            if(ntt.ge.90000)goto 202
            na2=na2+1
            read(s,8999)du,ia2(na2),du,aa2(na2),du,ra2(na2),du,
     &           chain2(na2),ir2(na2),du,xa2(na2),ya2(na2),za2(na2)
            ains2(na2)=du1
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.s(13:16).
     &           eq.'  CA')then
               i=i+1
               read(s,203)du,seqB(i),du,chB(i),nresB(i),du,
     &              xb(i),yb(i),zb(i)
               ins2(i)=du1
               do j=-1,20
                  if(seqB(i).eq.aa(j))then
                     seq1B(i)=slc(j)
                     goto 22
                  endif
               enddo
               seq1B(i)=slc(-1)
 22            continue
               if(i.ge.nmax)goto 202
            endif
         endif
      endif
      goto 201
 202  continue
 203  format(A17,A3,A1,A1,i4,A4,3F8.3)
 204  format(A100)
      close(10)
      nseqB=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

******************************************************************
*     pickup the aligned residues:
******************************************************************
      if(m_complex.eq.0)then    !monomer
         k=0
         do i=1,nseqA
            do j=1,nseqB
               if(nresA(i).eq.nresB(j))then
                  if(ins1(i).eq.ins2(j))then
                     k=k+1
                     iA(k)=i
                     iB(k)=j
                     goto 205
                  endif
               endif
            enddo
 205        continue
         enddo
      else                      !complex
         k=0
         do i=1,nseqA
            do j=1,nseqB
               if(nresA(i).eq.nresB(j).and.chA(i).eq.chB(j))then
                  if(ins1(i).eq.ins2(j))then !Residue_insertion_code
                     k=k+1
                     iA(k)=i
                     iB(k)=j
                     goto 206
                  endif
               endif
            enddo
 206        continue
         enddo
      endif
      n_ali=k                   !number of aligned residues
      if(n_ali.lt.1)then
         write(*,*)'There is no common residues in the input structures'
         goto 9999
      endif

******************************************************************
*     check the residue serial numbers ------------->
      if(m_complex.eq.0)then
         nA_repeat=0
         do i=1,nseqA
            do j=i+1,nseqA
               if(nresA(i).eq.nresA(j))then
                  if(ins1(i).eq.ins1(j))then
                     nA_repeat=nA_repeat+1
                  endif
               endif
            enddo
         enddo
         if(nA_repeat.gt.0)then
            write(*,380)nA_repeat
         endif
 380     format('Warning: TMscore calculation can be wrong because ',
     &        'there are ',I3,' residues with the serial number same ',
     &        'as others in first Chain. Please modify the PDB file ',
     &        'and rerun the program!!')
         nB_repeat=0
         do i=1,nseqB
            do j=i+1,nseqB
               if(nresB(i).eq.nresB(j))then
                  if(ins2(i).eq.ins2(j))then
                     nB_repeat=nB_repeat+1
                  endif
               endif
            enddo
         enddo
         if(nB_repeat.gt.0)then
            write(*,381)nB_repeat
         endif
 381     format('Warning: TMscore calculation can be wrong because ',
     &        'there are ',I3,' residues with the serial number same ',
     &        'as others in second Chain. Please modify the PDB file ',
     &        'and rerun the program!!')
      endif

************/////
*     parameters:
*****************
***   d0------------->
      if(nseqB.gt.15)then
         d0=1.24*(nseqB-15)**(1.0/3.0)-1.8
      else
         d0=0.5
      endif
      if(m_len.eq.1)then
         d0=1.24*(l0_fix-15)**(1.0/3.0)-1.8
      endif
      if(d0.lt.0.5)d0=0.5
      if(m_fix.eq.1)d0=d0_fix
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
***   iterative parameters ----->
      n_it=20                   !maximum number of iterations
      d_output=5                !for output alignment
      if(m_fix.eq.1)d_output=d0_fix
      n_init_max=6              !maximum number of L_init
      n_init=0
      L_ini_min=4
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

******************************************************************
*     find the maximum score starting from local structures superposition
******************************************************************
      score_max=-1              !TM-score
      score_maxsub_max=-1       !MaxSub-score
      score10_max=-1            !TM-score10
      n_GDT05_max=-1            !number of residues<0.5
      n_GDT1_max=-1             !number of residues<1
      n_GDT2_max=-1             !number of residues<2
      n_GDT4_max=-1             !number of residues<4
      n_GDT8_max=-1             !number of residues<8
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        do 300 iL=1,iL_max      !on aligned residues, [1,nseqA]
           LL=0
           ka=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              ka=ka+1
              k_ali(ka)=k
              LL=LL+1
           enddo
           if(i_init.eq.1)then  !global superposition
              call u3b(w,r_1,r_2,LL,2,rms,u,t,ier) !0:rmsd; 1:u,t; 2:rmsd,u,t
              armsd=dsqrt(rms/LL)
              rmsd_ali=armsd
           else
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun       !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
           if(score10_max.lt.score10)score10_max=score10
           if(score_maxsub_max.lt.score_maxsub)score_maxsub_max=
     &          score_maxsub
           if(n_GDT05_max.lt.n_GDT05)n_GDT05_max=n_GDT05
           if(n_GDT1_max.lt.n_GDT1)n_GDT1_max=n_GDT1
           if(n_GDT2_max.lt.n_GDT2)n_GDT2_max=n_GDT2
           if(n_GDT4_max.lt.n_GDT4)n_GDT4_max=n_GDT4
           if(n_GDT8_max.lt.n_GDT8)n_GDT8_max=n_GDT8
***   iteration for extending ---------------------------------->
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
              call score_fun    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(score10_max.lt.score10)score10_max=score10
              if(score_maxsub_max.lt.score_maxsub)score_maxsub_max
     &             =score_maxsub
              if(n_GDT05_max.lt.n_GDT05)n_GDT05_max=n_GDT05
              if(n_GDT1_max.lt.n_GDT1)n_GDT1_max=n_GDT1
              if(n_GDT2_max.lt.n_GDT2)n_GDT2_max=n_GDT2
              if(n_GDT4_max.lt.n_GDT4)n_GDT4_max=n_GDT4
              if(n_GDT8_max.lt.n_GDT8)n_GDT8_max=n_GDT8
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
 301       continue             !for iteration
 302       continue
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M

      ratio=1
      if(m_len.gt.0)then
         ratio=float(nseqB)/float(l0_fix)
      endif

******************************************************************
*     Output
******************************************************************
***   output TM-scale ---------------------------->
      write(*,*)
      write(*,*)'*****************************************************',
     &     '************************'
      write(*,*)'*                                 TM-SCORE           ',
     &     '                       *'
      write(*,*)'* A scoring function to assess the similarity of prot',
     &     'ein structures         *'
      write(*,*)'* Based on statistics:                               ',
     &     '                       *'
      write(*,*)'*       0.0 < TM-score < 0.17, random structural simi',
     &     'larity                 *'
      write(*,*)'*       0.5 < TM-score < 1.00, in about the same fold',
     &     '                       *'
      write(*,*)'* Reference: Yang Zhang and Jeffrey Skolnick, ',
     &     'Proteins 2004 57: 702-710     *'
      write(*,*)'* For comments, please email to: zhng@umich.edu      ',
     &     '                       *'
      write(*,*)'*****************************************************',
     &     '************************'
      write(*,*)
      write(*,501)pdb(1),nseqA
 501  format('Structure1: ',A10,'  Length= ',I4)
      if(m_len.eq.1)then
         write(*,411)pdb(2),nseqB
         write(*,412)l0_fix
      else
         write(*,502)pdb(2),nseqB
      endif
 411  format('Structure2: ',A10,'  Length= ',I4)
 412  format('TM-score is notmalized by ',I4)
 502  format('Structure2: ',A10,'  Length= ',I4,
     &     ' (by which all scores are normalized)')
      write(*,503)n_ali
 503  format('Number of residues in common= ',I4)
      write(*,513)rmsd_ali
 513  format('RMSD of  the common residues= ',F8.3)
      write(*,*)
      if(m_len.eq.1)then
         score_max=score_max*float(nseqB)/float(l0_fix)
      endif
      write(*,504)score_max,d0
 504  format('TM-score    = ',f6.4,'  (d0=',f5.2,')')
      write(*,505)score_maxsub_max*ratio
 505  format('MaxSub-score= ',f6.4,'  (d0= 3.50)')
      score_GDT=(n_GDT1_max+n_GDT2_max+n_GDT4_max+n_GDT8_max)
     &     /float(4*nseqB)
      write(*,506)score_GDT*ratio,n_GDT1_max/float(nseqB)*ratio,
     &     n_GDT2_max/float(nseqB)*ratio,n_GDT4_max/float(nseqB)*ratio,
     &     n_GDT8_max/float(nseqB)*ratio
 506  format('GDT-TS-score= ',f6.4,' %(d<1)=',f6.4,' %(d<2)=',f6.4,
     $     ' %(d<4)=',f6.4,' %(d<8)=',f6.4)
      score_GDT_HA=(n_GDT05_max+n_GDT1_max+n_GDT2_max+n_GDT4_max)
     &     /float(4*nseqB)
      write(*,507)score_GDT_HA*ratio,n_GDT05_max/float(nseqB)*ratio,
     &     n_GDT1_max/float(nseqB)*ratio,n_GDT2_max/float(nseqB)*ratio,
     &     n_GDT4_max/float(nseqB)*ratio
 507  format('GDT-HA-score= ',f6.4,' %(d<0.5)=',f6.4,' %(d<1)=',f6.4,
     $     ' %(d<2)=',f6.4,' %(d<4)=',f6.4)
      write(*,*)

***   recall and output the superposition of maxiumum TM-score:
      LL=0
      do i=1,ka0
         m=k_ali0(i)            !record of the best alignment
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,nseqA
         xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo

********* extract rotation matrix ------------>
      write(*,*)'-------- rotation matrix to rotate Chain-1 to ',
     &     'Chain-2 ------'
      write(*,*)'i          t(i)         u(i,1)         u(i,2) ',
     &     '        u(i,3)'
      do i=1,3
         write(*,304)i,t(i),u(i,1),u(i,2),u(i,3)
      enddo
c      do j=1,nseqA
c         xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
c         yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
c         zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
c         write(*,*)j,xt(j),yt(j),zt(j)
c      enddo
      write(*,*)
 304  format(I2,f18.10,f15.10,f15.10,f15.10)

********* rmsd in superposed regions --------------->
      d=d_output                !for output
      call score_fun()          !give i_ali(i), score_max=score now

***   output rotated chain1 + chain2----->
      if(m_out.ne.1)goto 999
ccc   output CA-trace superposition------>
      OPEN(unit=7,file=outname,status='unknown') !pdb1.aln + pdb2.aln
 900  format(A)
 901  format('select atomno= ',I5)
 902  format('select atomno= ',I5)
      write(7,900)'load inline'
      write(7,900)'select atomno<50001'
      write(7,900)'wireframe .45'
      write(7,900)'select none'
      write(7,900)'select atomno>50000'
      write(7,900)'wireframe .15'
      write(7,900)'color white' !all above atoms are white
      do i=1,n_cut
         write(7,901)iA(i_ali(i))
         write(7,900)'color red'
         write(7,902)50000+iB(i_ali(i))
         write(7,900)'color red'
      enddo
      write(7,900)'select all'
      write(7,900)'exit'
      write(7,514)rmsd_ali
 514  format('REMARK  RMSD of the common residues=',F8.3)
      write(7,515)score_max,d0
 515  format('REMARK  TM-score=',f6.4,' (d0=',f5.2,')')
      do i=1,nseqA
         write(7,1236)i,seqA(i),chA(i),nresA(i),ins1(i),
     &        xt(i),yt(i),zt(i)
      enddo
      write(7,1238)
      do i=2,nseqA
         write(7,1239)i-1,i
      enddo
      do i=1,nseqB
         write(7,1236)50000+i,seqB(i),chB(i),nresB(i),ins2(i),
     &        xb(i),yb(i),zb(i)
      enddo
      write(7,1238)
      do i=2,nseqB
         write(7,1239)50000+i-1,50000+i
      enddo
      close(7)
 1236 format('ATOM  ',i5,'  CA  ',A3,' ',A1,I4,A1,3X,3F8.3)
 1238 format('TER')
 1239 format('CONECT',I5,I5)
ccc   output all-atom superposition------>
c      outname=trim(outname)//'_atm'
      outname=outname(1:len_trim(outname))//'_atm'
      OPEN(unit=7,file=outname,status='unknown') !pdb1.aln + pdb2.aln
      write(7,900)'load inline'
      write(7,900)'select atomno<50001'
      write(7,900)'color blue'
      write(7,900)'select atomno>50000'
      write(7,900)'color red'
      write(7,900)'select all'
      write(7,900)'cartoon'
      write(7,900)'exit'
      write(7,514)rmsd_ali
      write(7,515)score_max,d0
***   chain1:
      do i=1,na1
         ax=t(1)+u(1,1)*xa1(i)+u(1,2)*ya1(i)+u(1,3)*za1(i)
         ay=t(2)+u(2,1)*xa1(i)+u(2,2)*ya1(i)+u(2,3)*za1(i)
         az=t(3)+u(3,1)*xa1(i)+u(3,2)*ya1(i)+u(3,3)*za1(i)
         write(7,8888)i,aa1(i),ra1(i),chain1(i),ir1(i),ains1(i),ax,ay,az
      enddo
      write(7,1238)             !TER
***   chain2:
      do i=1,na2
         write(7,8888)50000+i,aa2(i),ra2(i),chain2(i),ir2(i),ains2(i),
     &        xa2(i),ya2(i),za2(i)
      enddo
      write(7,1238)             !TER
      close(7)
 8888 format('ATOM  ',I5,1x,A4,1x,A3,' ',A1,I4,A1,3x,3F8.3)
 999  continue

***   record aligned residues by i=[1,nseqA], for sequenceM()------------>
      do i=1,nseqA
         iq(i)=0
      enddo
      do i=1,n_cut
         j=iA(i_ali(i))         ![1,nseqA]
         k=iB(i_ali(i))         ![1,nseqB]
         dis=sqrt((xt(j)-xb(k))**2+(yt(j)-yb(k))**2+(zt(j)-zb(k))**2)
c         write(*,*)i,j,k,dis,d_output,'1--'
         if(dis.lt.d_output)then
            iq(j)=1
         endif
      enddo
*******************************************************************
***   output aligned sequences
      if(m_complex.eq.0)then
         k=0
         i=1
         j=1
 800     continue
         if(i.gt.nseqA.and.j.gt.nseqB)goto 802
         if(i.gt.nseqA.and.j.le.nseqB)then
            k=k+1
            sequenceA(k)='-'
            sequenceB(k)=seq1B(j)
            sequenceM(k)=' '
            j=j+1
            goto 800
         endif
         if(i.le.nseqA.and.j.gt.nseqB)then
            k=k+1
            sequenceA(k)=seq1A(i)
            sequenceB(k)='-'
            sequenceM(k)=' '
            i=i+1
            goto 800
         endif
         if(nresA(i).eq.nresB(j))then
            k=k+1
            sequenceA(k)=seq1A(i)
            sequenceB(k)=seq1B(j)
            if(iq(i).eq.1)then
               sequenceM(k)=':'
            else
               sequenceM(k)=' '
            endif
            i=i+1
            j=j+1
            goto 800
         elseif(nresA(i).lt.nresB(j))then
            k=k+1
            sequenceA(k)=seq1A(i)
            sequenceB(k)='-'
            sequenceM(k)=' '
            i=i+1
            goto 800
         elseif(nresB(j).lt.nresA(i))then
            k=k+1
            sequenceA(k)='-'
            sequenceB(k)=seq1B(j)
            sequenceM(k)=' '
            j=j+1
            goto 800
         endif
 802     continue
      else
         k=0
         i=1                    !
         j=1
 803     continue
         if(i.gt.nseqA.and.j.gt.nseqB)goto 804
         if(i.gt.nseqA.and.j.le.nseqB)then
            k=k+1
            sequenceA(k)='-'
            sequenceB(k)=seq1B(j)
            sequenceM(k)=' '
            j=j+1
            goto 803
         endif
         if(i.le.nseqA.and.j.gt.nseqB)then
            k=k+1
            sequenceA(k)=seq1A(i)
            sequenceB(k)='-'
            sequenceM(k)=' '
            i=i+1
            goto 803
         endif
         if(chA(i).ne.chB(j))then
            kk=0                !check if chA(i) match later chains
            do i2=j,nseqB
               if(chA(i).eq.chB(i2))then
                  kk=kk+1
               endif
            enddo
            if(kk.eq.0)then     !chA(i) is not matched
               k=k+1
               sequenceA(k)=seq1A(i)
               sequenceB(k)='-'
               sequenceM(k)=' '
               i=i+1
               goto 803
            endif

            kk=0                !check if chB(j) match later chains
            do i1=i,nseqA
               if(chA(i1).eq.chB(j))then
                  kk=kk+1
               endif
            enddo
            if(kk.eq.0)then     !chB(j) is not matched
               k=k+1
               sequenceA(k)='-'
               sequenceB(k)=seq1B(j)
               sequenceM(k)=' '
               j=j+1
               goto 803
            endif
            write(*,*)'Chains in complexes have crossed order, ',
     &           'therefore, no alignment is output'
            stop
         else
            if(nresA(i).eq.nresB(j))then
               k=k+1
               sequenceA(k)=seq1A(i)
               sequenceB(k)=seq1B(j)
               if(iq(i).eq.1)then
                  sequenceM(k)=':'
               else
                  sequenceM(k)=' '
               endif
               i=i+1
               j=j+1
               goto 803
            elseif(nresA(i).lt.nresB(j))then
               k=k+1
               sequenceA(k)=seq1A(i)
               sequenceB(k)='-'
               sequenceM(k)=' '
               i=i+1
               goto 803
            elseif(nresB(j).lt.nresA(i))then
               k=k+1
               sequenceA(k)='-'
               sequenceB(k)=seq1B(j)
               sequenceM(k)=' '
               j=j+1
               goto 803
            endif
         endif
 804     continue
      endif

ccc   RMSD (d<5.0)-------->
      LL=0
      do i=1,n_cut
         m=i_ali(i)             ![1,nseqA]
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,0,rms,u,t,ier)
      armsd=dsqrt(rms/LL)
      rmsd=armsd

      write(*,600)d_output,n_cut,rmsd
 600  format('Superposition in the TM-score: Length(d<',f3.1,
     $     ')=',i3,'  RMSD=',f6.2)
      write(*,603)d_output
 603  format('(":" denotes the residue pairs of distance < ',f3.1,
     &     ' Angstrom)')
      write(*,601)(sequenceA(i),i=1,k)
      write(*,601)(sequenceM(i),i=1,k)
      write(*,601)(sequenceB(i),i=1,k)
      write(*,602)(mod(i,10),i=1,k)
 601  format(2000A1)
 602  format(2000I1)
      write(*,*)

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 9999 END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine score_fun
      PARAMETER(nmax=5000)

      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0,d0_fix
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score,score_maxsub,score_fix,score10
      common/GDT/n_GDT05,n_GDT1,n_GDT2,n_GDT4,n_GDT8
      double precision score,score_max,score_fix,score_fix_max
      double precision score_maxsub,score10

      d_tmp=d
 21   n_cut=0                   !number of residue-pairs dis<d, for iteration
      n_GDT05=0                 !for GDT-score, # of dis<0.5
      n_GDT1=0                  !for GDT-score, # of dis<1
      n_GDT2=0                  !for GDT-score, # of dis<2
      n_GDT4=0                  !for GDT-score, # of dis<4
      n_GDT8=0                  !for GDT-score, # of dis<8
      score_maxsub_sum=0        !Maxsub-score
      score_sum=0               !TMscore
      score_sum10=0             !TMscore10
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis=sqrt((xt(i)-xb(j))**2+(yt(i)-yb(j))**2+(zt(i)-zb(j))**2)
***   for iteration:
         if(dis.lt.d_tmp)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
***   for GDT-score:
         if(dis.le.8)then
            n_GDT8=n_GDT8+1
            if(dis.le.4)then
               n_GDT4=n_GDT4+1
               if(dis.le.2)then
                  n_GDT2=n_GDT2+1
                  if(dis.le.1)then
                     n_GDT1=n_GDT1+1
                     if(dis.le.0.5)then
                        n_GDT05=n_GDT05+1
                     endif
                  endif
               endif
            endif
         endif
***   for MAXsub-score:
         if(dis.lt.3.5)then
            score_maxsub_sum=score_maxsub_sum+1/(1+(dis/3.5)**2)
         endif
***   for TM-score:
         score_sum=score_sum+1/(1+(dis/d0)**2)
***   for TM-score10:
         if(dis.lt.10)then
            score_sum10=score_sum10+1/(1+(dis/d0)**2)
         endif
      enddo
      if(n_cut.lt.3.and.n_ali.gt.3)then
         d_tmp=d_tmp+.5
         goto 21
      endif
      score_maxsub=score_maxsub_sum/float(nseqB) !MAXsub-score
      score=score_sum/float(nseqB) !TM-score
      score10=score_sum10/float(nseqB) !TM-score10

      return
      end

cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
c  w    - w(m) is weight for atom pair  c m                    (given)
c  x    - x(i,m) are coordinates of atom c m in set x          (given)
c  y    - y(i,m) are coordinates of atom c m in set y          (given)
c  n    - n is number of atom pairs                            (given)
c  mode  - 0:calculate rms     only                            (given,short)
c          1:calculate     u,t only                            (given,medium)
c          2:calculate rms,u,t                                 (given,longer)
c  rms   - sum of w*(ux+t-y)**2 over all atom pairs            (result)
c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
c  t    - t(i)   is translation vector for best superposition  (result)
c  ier  - 0: a unique optimal superposition has been determined(result)
c       -1: superposition is not unique but optimal
c       -2: no result obtained because of negative weights w
c           or all weights equal to zero.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
      double precision w(*), x(3,*), y(3,*)
      integer n, mode

      double precision rms, u(3,3), t(3)
      integer ier

      integer i, j, k, l, m1, m
      integer ip(9), ip2312(4)
      double precision r(3,3), xc(3), yc(3), wc
      double precision a(3,3), b(3,3), e(3), rr(6), ss(6)
      double precision e0, d, spur, det, cof, h, g
      double precision cth, sth, sqrth, p, sigma
      double precision c1x, c1y, c1z, c2x, c2y, c2z
      double precision s1x, s1y, s1z, s2x, s2y, s2z
      double precision sxx, sxy, sxz, syx, syy, syz, szx, szy, szz

      double precision sqrt3, tol, zero

      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /

      wc  = zero
      rms = zero
      e0  = zero
      s1x = zero
      s1y = zero
      s1z = zero
      s2x = zero
      s2y = zero
      s2z = zero
      sxx = zero
      sxy = zero
      sxz = zero
      syx = zero
      syy = zero
      syz = zero
      szx = zero
      szy = zero
      szz = zero

      do i=1, 3
         xc(i) = zero
         yc(i) = zero
         t(i) = zero
         do j=1, 3
            r(i,j) = zero
            u(i,j) = zero
            a(i,j) = zero
            if( i .eq. j ) then
               u(i,j) = 1.0
               a(i,j) = 1.0
            end if
         end do
      end do

      ier = -1
      if( n .lt. 1 ) return
      ier = -2

      do m=1, n
         c1x=x(1, m)
         c1y=x(2, m)
         c1z=x(3, m)

         c2x=y(1, m)
         c2y=y(2, m)
         c2z=y(3, m)

         s1x = s1x + c1x
         s1y = s1y + c1y;
         s1z = s1z + c1z;

         s2x = s2x + c2x;
         s2y = s2y + c2y;
         s2z = s2z + c2z;

         sxx = sxx + c1x*c2x;
         sxy = sxy + c1x*c2y;
         sxz = sxz + c1x*c2z;

         syx = syx + c1y*c2x;
         syy = syy + c1y*c2y;
         syz = syz + c1y*c2z;

         szx = szx + c1z*c2x;
         szy = szy + c1z*c2y;
         szz = szz + c1z*c2z;
      end do

      xc(1) = s1x/n;
      xc(2) = s1y/n;
      xc(3) = s1z/n;

      yc(1) = s2x/n;
      yc(2) = s2y/n;
      yc(3) = s2z/n;
      if(mode.eq.2.or.mode.eq.0) then ! need rmsd
         do m=1, n
            do i=1, 3
               e0 = e0+ (x(i, m)-xc(i))**2 + (y(i, m)-yc(i))**2
            end do
         end do
      endif

      r(1, 1) = sxx-s1x*s2x/n;
      r(2, 1) = sxy-s1x*s2y/n;
      r(3, 1) = sxz-s1x*s2z/n;
      r(1, 2) = syx-s1y*s2x/n;
      r(2, 2) = syy-s1y*s2y/n;
      r(3, 2) = syz-s1y*s2z/n;
      r(1, 3) = szx-s1z*s2x/n;
      r(2, 3) = szy-s1z*s2y/n;
      r(3, 3) = szz-s1z*s2z/n;

      det = r(1,1) * ( (r(2,2)*r(3,3)) - (r(2,3)*r(3,2)) )
     &     - r(1,2) * ( (r(2,1)*r(3,3)) - (r(2,3)*r(3,1)) )
     &     + r(1,3) * ( (r(2,1)*r(3,2)) - (r(2,2)*r(3,1)) )

      sigma = det

      m = 0
      do j=1, 3
         do i=1, j
            m = m+1
            rr(m) = r(1,i)*r(1,j) + r(2,i)*r(2,j) + r(3,i)*r(3,j)
         end do
      end do

      spur = (rr(1)+rr(3)+rr(6)) / 3.0
      cof = (((((rr(3)*rr(6) - rr(5)*rr(5)) + rr(1)*rr(6))
     &     - rr(4)*rr(4)) + rr(1)*rr(3)) - rr(2)*rr(2)) / 3.0
      det = det*det

      do i=1, 3
         e(i) = spur
      end do
      if( spur .le. zero ) goto 40
      d = spur*spur
      h = d - cof
      g = (spur*cof - det)/2.0 - spur*h
      if( h .le. zero ) then
         if( mode .eq. 0 ) then
            goto 50
         else
            goto 30
         end if
      end if
      sqrth = dsqrt(h)
      d = h*h*h - g*g
      if( d .lt. zero ) d = zero
      d = datan2( dsqrt(d), -g ) / 3.0
      cth = sqrth * dcos(d)
      sth = sqrth*sqrt3*dsin(d)
      e(1) = (spur + cth) + cth
      e(2) = (spur - cth) + sth
      e(3) = (spur - cth) - sth

      if( mode .eq. 0 ) then
         goto 50
      end if

      do l=1, 3, 2
         d = e(l)
         ss(1) = (d-rr(3)) * (d-rr(6))  - rr(5)*rr(5)
         ss(2) = (d-rr(6)) * rr(2)      + rr(4)*rr(5)
         ss(3) = (d-rr(1)) * (d-rr(6))  - rr(4)*rr(4)
         ss(4) = (d-rr(3)) * rr(4)      + rr(2)*rr(5)
         ss(5) = (d-rr(1)) * rr(5)      + rr(2)*rr(4)
         ss(6) = (d-rr(1)) * (d-rr(3))  - rr(2)*rr(2)

         if( dabs(ss(1)) .ge. dabs(ss(3)) ) then
            j=1
            if( dabs(ss(1)) .lt. dabs(ss(6)) ) j = 3
         else if( dabs(ss(3)) .ge. dabs(ss(6)) ) then
            j = 2
         else
            j = 3
         end if

         d = zero
         j = 3 * (j - 1)

         do i=1, 3
            k = ip(i+j)
            a(i,l) = ss(k)
            d = d + ss(k)*ss(k)
         end do
         if( d .gt. zero ) d = 1.0 / dsqrt(d)
         do i=1, 3
            a(i,l) = a(i,l) * d
         end do
      end do

      d = a(1,1)*a(1,3) + a(2,1)*a(2,3) + a(3,1)*a(3,3)
      if ((e(1) - e(2)) .gt. (e(2) - e(3))) then
         m1 = 3
         m = 1
      else
         m1 = 1
         m = 3
      endif

      p = zero
      do i=1, 3
         a(i,m1) = a(i,m1) - d*a(i,m)
         p = p + a(i,m1)**2
      end do
      if( p .le. tol ) then
         p = 1.0
         do 21 i=1, 3
            if (p .lt. dabs(a(i,m))) goto 21
            p = dabs( a(i,m) )
            j = i
 21      continue
         k = ip2312(j)
         l = ip2312(j+1)
         p = dsqrt( a(k,m)**2 + a(l,m)**2 )
         if( p .le. tol ) goto 40
         a(j,m1) = zero
         a(k,m1) = -a(l,m)/p
         a(l,m1) =  a(k,m)/p
      else
         p = 1.0 / dsqrt(p)
         do i=1, 3
            a(i,m1) = a(i,m1)*p
         end do
      end if

      a(1,2) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
      a(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
      a(3,2) = a(1,3)*a(2,1) - a(1,1)*a(2,3)

 30   do l=1, 2
         d = zero
         do i=1, 3
            b(i,l) = r(i,1)*a(1,l) + r(i,2)*a(2,l) + r(i,3)*a(3,l)
            d = d + b(i,l)**2
         end do
         if( d .gt. zero ) d = 1.0 / dsqrt(d)
         do i=1, 3
            b(i,l) = b(i,l)*d
         end do
      end do
      d = b(1,1)*b(1,2) + b(2,1)*b(2,2) + b(3,1)*b(3,2)
      p = zero

      do i=1, 3
         b(i,2) = b(i,2) - d*b(i,1)
         p = p + b(i,2)**2
      end do
      if( p .le. tol ) then
         p = 1.0
         do 22 i=1, 3
            if(p.lt.dabs(b(i,1)))goto 22
            p = dabs( b(i,1) )
            j = i
 22      continue
         k = ip2312(j)
         l = ip2312(j+1)
         p = dsqrt( b(k,1)**2 + b(l,1)**2 )
         if( p .le. tol ) goto 40
         b(j,2) = zero
         b(k,2) = -b(l,1)/p
         b(l,2) =  b(k,1)/p
      else
         p = 1.0 / dsqrt(p)
         do i=1, 3
            b(i,2) = b(i,2)*p
         end do
      end if

      b(1,3) = b(2,1)*b(3,2) - b(2,2)*b(3,1)
      b(2,3) = b(3,1)*b(1,2) - b(3,2)*b(1,1)
      b(3,3) = b(1,1)*b(2,2) - b(1,2)*b(2,1)

      do i=1, 3
         do j=1, 3
            u(i,j) = b(i,1)*a(j,1) + b(i,2)*a(j,2) + b(i,3)*a(j,3)
         end do
      end do

 40   do i=1, 3
         t(i) = ((yc(i) - u(i,1)*xc(1)) - u(i,2)*xc(2)) - u(i,3)*xc(3)
      end do
 50   do i=1, 3
         if( e(i) .lt. zero ) e(i) = zero
         e(i) = dsqrt( e(i) )
      end do

      ier = 0
      if( e(2) .le. (e(1) * 1.0d-05) ) ier = -1

      d = e(3)
      if( sigma .lt. 0.0 ) then
         d = - d
         if( (e(2) - e(3)) .le. (e(1) * 1.0d-05) ) ier = -1
      end if
      d = (d + e(2)) + e(1)

      if(mode .eq. 2.or.mode.eq.0) then ! need rmsd
         rms = (e0 - d) - d
         if( rms .lt. 0.0 ) rms = 0.0
      endif

      return
      end


