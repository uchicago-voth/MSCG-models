!
!Copyright 2015 University of Chicago
!
!Licensed under the Apache License, Version 2.0 (the "License");
!you may not use this file except in compliance with the License.
!You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
!Unless required by applicable law or agreed to in writing, software
!distributed under the License is distributed on an "AS IS" BASIS,
!WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!See the License for the specific language governing permissions and
!limitations under the License.
!
! Provides definitions of linked lists.
!

module obj_ll

        use env_kindtypes, only: si

        private

        public          :: i_sp_dLL

        type :: i_sp_dLL_link
                private

                type    (i_sp_dLL_link),pointer:: next   => null()
                type    (i_sp_dLL_link),pointer:: prev   => null()
                integer (si)                   :: val
        end type i_sp_dLL_link

        type :: i_sp_dLL
                private

                type    (i_sp_dLL_link),pointer :: head_link    => null()
                type    (i_sp_dLL_link),pointer :: current_link => null()
                type    (i_sp_dLL_link),pointer :: tail_link    => null()
                integer (si)                    :: length = 0

                contains
                        procedure,public   :: initialize=> i_sp_dLL_init
                        procedure,public   :: prepend   => i_sp_dLL_prepend 
                        procedure,public   :: append    => i_sp_dLL_append 
                        procedure,public   :: insert    => i_sp_dLL_insert
                        procedure,public   :: delete    => i_sp_dLL_delete
                        procedure,public   :: destroy   => i_sp_dLL_destroy

                        procedure,public   :: getLength => i_sp_dLL_getLength
                        procedure,public   :: getAll    => i_sp_dLL_getAll
                        procedure,public   :: sgetAll   => i_sp_dLL_sgetAll
                        procedure,public   :: getCurrent=> i_sp_dLL_getCurrent

                        procedure,public   :: next      => i_sp_dLL_incrCurr
                        procedure,public   :: prev      => i_sp_dLL_decrCurr
                        procedure,public   :: goto      => i_sp_dLL_goto

                        procedure,public   :: presentQuery => i_sp_dLL_presentQuery

                        final              :: i_sp_dLL_final
        end type i_sp_dLL

        contains
                pure subroutine i_sp_dLL_init(self,toConvert)

                        implicit none

                        class(i_sp_dLL),intent(inout)           :: self
                        integer (si),   intent(in   ),optional  :: toConvert(:)

                        integer (si)                  :: iter

                        if (self%length /= 0) then
                                call self%destroy()
                        endif

                        if (present(toConvert)) then
                                do iter=1,size(toConvert)
                                        call self%append(toConvert(iter))
                                enddo
                        endif

                endsubroutine i_sp_dLL_init

                elemental subroutine i_sp_dLL_destroy(self)

                        implicit none

                        class(i_sp_dLL),intent(inout) :: self

                        type(i_sp_dLL_link),pointer   :: placeholder
                        integer (si)                  :: iter

                        if (self%length==0) then
                                !nothing to be done
                                return
                        endif

                        call self%goto(1)
                        do iter=1,(self%length-1)
                                placeholder => self%current_link%next
                                deallocate(self%current_link)
                                self%current_link => placeholder
                        enddo
                        !last link
                        deallocate(self%current_link)

                        nullify(self%head_link)
                        nullify(self%tail_link)

                        self%length = 0

                endsubroutine i_sp_dLL_destroy

                elemental subroutine i_sp_dLL_final(self)

                        implicit none

                        type(i_sp_dLL),intent(inout) :: self

                        type(i_sp_dLL_link),pointer   :: placeholder
                        integer (si)                  :: iter

                        if (self%length==0) then
                                !nothing to be done
                                return
                        endif

                        call self%goto(1)
                        do iter=1,(self%length-1)
                                placeholder => self%current_link%next
                                deallocate(self%current_link)
                                self%current_link => placeholder
                        enddo
                        !last link
                        deallocate(self%current_link)

                        nullify(self%head_link)
                        nullify(self%tail_link)

                        self%length = 0

                endsubroutine i_sp_dLL_final

                pure subroutine i_sp_dLL_append(self,valToAdd)

                        implicit none

                        class(i_sp_dLL),intent(inout) :: self
                        integer (si)   ,intent(in   ) :: valToAdd

                        if (self%length==0) then
                                !we're an empty list
                                allocate(self%head_link)
                                self%tail_link     => self%head_link
                                self%tail_link%val =  valToAdd

                                !init current place
                                self%current_link => self%head_link
                        else
                                !we're not an empty list
                                allocate(self%tail_link%next)
                                self%tail_link%next%prev => self%tail_link
                                self%tail_link => self%tail_link%next
                                self%tail_link%val = valToAdd
                        endif

                        self%length = self%length + 1
                endsubroutine i_sp_dLL_append

                pure subroutine i_sp_dLL_prepend(self,valToAdd)

                        implicit none

                        class(i_sp_dLL),intent(inout) :: self
                        integer (si)   ,intent(in   ) :: valToAdd

                        if (self%length==0) then
                                !we're an empty list
                                allocate(self%head_link)
                                self%tail_link     => self%head_link
                                self%head_link%val =  valToAdd

                                !init current place
                                self%current_link => self%head_link
                        else
                                !we're not an empty list
                                allocate(self%head_link%prev)
                                self%head_link%prev%next => self%head_link
                                self%head_link => self%head_link%prev
                                self%head_link%val = valToAdd
                        endif

                        self%length = self%length + 1
                endsubroutine i_sp_dLL_prepend

                pure subroutine i_sp_dLL_insert(self,valToAdd,after)
                !there's a non polymorphic interval varaible here. 
                !in the future this may cause issues...?

                        implicit none

                        class(i_sp_dLL),intent(inout)           :: self
                        integer (si)   ,intent(in   )           :: valToAdd
                        logical,        intent(in   ),optional  :: after

                        type(i_sp_dLL_link),pointer   :: placeholder
                        logical                       :: after_

                        !default arg guard
                        if (present(after)) then
                                after_ = after
                        else
                                after_ = .false.
                        endif

                        if (self%length==0) then
                                !we're an empty list
                                allocate(self%head_link)
                                self%tail_link     => self%head_link
                                self%head_link%val =  valToAdd

                                !init current place
                                self%current_link => self%head_link
                        elseif (after_) then
                                if (associated(self%current_link,self%tail_link)) then
                                        call self%append(valToAdd)
                                        return
                                endif
                                !not at end
                                placeholder => self%current_link%next

                                nullify(self%current_link%next)

                                !this now points to the new node.
                                allocate(self%current_link%next)

                                self%current_link%next%prev => self%current_link
                                placeholder%prev            => self%current_link%next
                                self%current_link%next%next => placeholder

                                self%current_link%next%val = valToAdd
                        else
                                if (associated(self%current_link,self%head_link)) then
                                        call self%prepend(valToAdd)
                                        return
                                endif
                                !not at beginning
                                placeholder => self%current_link%prev

                                nullify(self%current_link%prev)

                                !this now points to the new node.
                                allocate(self%current_link%prev)

                                self%current_link%prev%next => self%current_link
                                placeholder%next            => self%current_link%prev
                                self%current_link%prev%prev => placeholder

                                self%current_link%prev%val = valToAdd
                        endif

                        self%length = self%length + 1
                endsubroutine i_sp_dLL_insert

                pure subroutine i_sp_dLL_delete(self)
                !there's a non polymorphic interval varaible here. 
                !in the future this may cause issues...?

                        implicit none

                        class(i_sp_dLL),intent(inout) :: self

                        type(i_sp_dLL_link),pointer   :: placeholder

                        if (self%length==0) then
                                !we're an empty list
                                return
                        else if (self%length==1) then
                                deallocate(self%head_link)
                                self%current_link => null()
                                self%tail_link => null()
                        else
                                !we're not an empty list
                                placeholder => self%current_link

                                if (.not. associated(self%current_link%prev)) then
                                        !at the head node
                                        self%current_link%next%prev => null()
                                        self%head_link => self%current_link%next
                                        self%current_link => self%head_link
                                elseif (.not. associated(self%current_link%next)) then
                                        !at the tail node
                                        self%current_link%prev%next => null()
                                        self%tail_link => self%current_link%prev
                                        !our norm is to point to the *next* node
                                        self%current_link => self%head_link
                                else
                                        self%current_link%prev%next => self%current_link%next
                                        self%current_link%next%prev => self%current_link%prev
                                        self%current_link           => self%current_link%next
                                endif

                                deallocate(placeholder)
                        endif

                        self%length = self%length - 1
                endsubroutine i_sp_dLL_delete

                pure subroutine i_sp_dll_incrCurr(self)

                        implicit none

                        class(i_sp_dll),intent(inout) :: self

                        if (self%length==0) then
                                !we're an empty list
                                return
                        elseif (associated(self%current_link,self%tail_link)) then
                                !wrap around to beginning
                                self%current_link => self%head_link
                        else
                                self%current_link => self%current_link%next
                        endif
                endsubroutine i_sp_dll_incrCurr

                pure subroutine i_sp_dll_decrCurr(self)

                        implicit none

                        class(i_sp_dll),intent(inout) :: self

                        if (self%length==0) then
                                !we're an empty list
                                return
                        elseif (associated(self%current_link,self%head_link)) then
                                !wrap around to beginning
                                self%current_link => self%tail_link
                        else
                                self%current_link => self%current_link%prev
                        endif
                endsubroutine i_sp_dll_decrCurr

                elemental subroutine i_sp_dll_goto(self,which)

                        implicit none

                        class(i_sp_dll),intent(inout) :: self
                        integer (si)   ,intent(in   ) :: which

                        integer (si)                  :: iter

                        if (self%length==0) then
                                !we're an empty list
                                return
                        elseif (which==1) then
                                !goto first element
                                self%current_link => self%head_link
                                return
                        elseif (which>=self%length) then
                                !goto last element, even if overshot
                                self%current_link => self%tail_link
                                return
                        elseif (which>=self%length/2) then
                                !goto last element, travel
                                self%current_link => self%tail_link
                                do iter=1,(self%length - which)
                                        self%current_link => self%current_link%prev
                                enddo
                        else
                                !goto first element, travel
                                self%current_link => self%head_link
                                do iter=1,(which-1)
                                        self%current_link => self%current_link%next
                                enddo
                        endif
                endsubroutine i_sp_dll_goto

                elemental function i_sp_dll_getCurrent(self)
                        !this can segfault. if there are no elements

                        implicit none

                        class(i_sp_dll),intent(in   ) :: self
                        integer (si)                  :: i_sp_dLL_getCurrent

                        i_sp_dLL_getCurrent = self%current_link%val

                endfunction i_sp_dll_getCurrent

                elemental function i_sp_dll_getLength(self)

                        implicit none

                        class(i_sp_dll),intent(in   ) :: self
                        integer (si)                  :: i_sp_dLL_getLength

                        if (self%length==0) then
                                !we're an empty list
                                i_sp_dLL_getLength = 0
                        else
                                i_sp_dLL_getLength = self%length
                        endif
                endfunction i_sp_dll_getLength

                pure subroutine i_sp_dll_sgetAll(self,toReturn)

                        implicit none

                        class(i_sp_dll),intent(inout)             :: self
                        integer (si),   intent(inout),allocatable :: toReturn(:)

                        type(i_sp_dll_link),pointer   :: iterator
                        integer (si)                  :: iter

                        if (self%length==0) then
                                !we're an empty list

                                !careful allocation
                                if (allocated(toReturn)) then
                                        if (.not. (size(toReturn,1) .eq. 0)) then
                                                deallocate(toReturn)
                                                allocate(toReturn(0))
                                        endif
                                else
                                        allocate(toReturn(0))
                                endif
                                return
                        else
                                !careful allocation
                                if (allocated(toReturn)) then
                                        if (.not. (size(toReturn,1) .eq. self%length)) then
                                                deallocate(toReturn)
                                                allocate(toReturn(self%length))
                                        endif
                                else
                                        allocate(toReturn(self%length))
                                endif
                                iterator => self%head_link
                                do iter=1,self%length
                                        toReturn(iter) = iterator%val
                                        iterator => iterator%next
                                enddo
                        endif
                endsubroutine i_sp_dll_sgetAll

                function i_sp_dll_getAll(self) result(toReturn)

                        implicit none

                        class(i_sp_dll),intent(inout)             :: self
                        integer (si),                 allocatable :: toReturn(:)

                        call i_sp_dll_sgetAll(self,toReturn)

                endfunction i_sp_dll_getAll

                function i_sp_dll_presentQuery(self,query) result(toReturn)

                        implicit none

                        class(i_sp_dll),intent(inout)   :: self
                        integer (si),   intent(in   )   :: query

                        logical                         :: toReturn

                        type(i_sp_dll_link),pointer     :: iterator
                        integer (si)                    :: iter

                        toReturn = .false.

                        if (self%length==0) then
                                !we're an empty list
                                return
                        else
                                iterator => self%head_link
                                find: do iter=1,self%getLength()
                                        if (iterator%val == query) then
                                                toReturn = .true.
                                                exit find
                                        else
                                                iterator => iterator%next
                                        endif
                                enddo find
                                return
                        endif
                endfunction i_sp_dll_presentQuery
endmodule obj_ll

!program test
!
!        use env_kindtypes,      only: si
!        use obj_ll,             only: i_sp_dLL
!
!        type(i_sp_dLL)     :: list
!
!        integer (si)       :: iter
!        integer (si)       :: length = 10
!
!        call list%initialize()
!        do iter=1,length
!                call list%insert(iter,.false.)
!                print*, list%getCurrent(), list%getLength()
!        enddo
!        print*, list%getAll()
!        call list%destroy()
!
!        call list%initialize([ 0 , 1 , 2 , 3 ])
!        print*, list%getAll()
!        print*, list%getCurrent()
!        call    list%next()
!        print*, list%getCurrent()
!        call    list%prev()
!        print*, list%getCurrent()
!        call    list%prev()
!        print*, list%getCurrent()
!        call    list%prev()
!        print*, list%getCurrent()
!        call    list%goto(2)
!        print*, list%getCurrent()
!        call    list%delete()
!        print*, list%getCurrent()
!        print*, list%getAll()
!        call    list%prepend(9)
!        call    list%append(8)
!        print*, list%getAll()
!        print*, list%presentQuery(9)
!        print*, list%presentQuery(10)
!
!endprogram
!
