
.. _governance:

PypeIt Governance Structure
===========================

PypeIt follows a structured do-ocracy governance model, where valuable and
sustained contributions to the project lead to more influence over its
development direction.

.. note::

    We are currently in a transition to a new governance model that evolves the
    current structure into one that is intended to be more community focused,
    but remains primarily a do-ocracy.  We expect the new governance model will
    be instituted by the end of 2025.  The following is **a working draft** of
    this governance structure.

==========   ===========  ============= ====================
Date         Status       Authors       Comments
==========   ===========  ============= ====================
2025-07-29   DRAFT        Westfall      First draft
==========   ===========  ============= ====================

.. contents:: Sections
    :depth: 2
    :local:

.. _governance-structure:

Project Committees and Team Roles
---------------------------------

The PypeIt governance model establishes the following committees and team roles,
as defined by this section.  The process of appointing each role is discussed
:ref:`below<governance-appointments>`.

- :ref:`Coordination Committee<governance-coordination>`: Provides overall
  project management and coordinates development efforts by setting goals and
  milestones that incorporate community input.

- :ref:`Advisory Council<governance-advisory>`: A self-governing body that, in
  part, provides outward-focused reviews of PypeIt's impact and success, as
  needed to ensure PypeIt is able to meet the needs of the larger astronomy
  community.  See their full charter :ref:`here<advisory-council-charter>`.

- :ref:`PypeIt Maintainers<governance-maintain>`: Individuals with write access
  to the main repository that supervise and contribute to code development.

- :ref:`PypeIt Voting Members<governance-voting>`: Individuals that vote on
  project-related decisions on behalf of the community, keeping in mind the best
  interests of both the broader community and the PypeIt project.
  
- :ref:`Other Roles<governance-other>`

.. _governance-coordination:

PypeIt Coordination Committee
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The PypeIt Coordination Committee (PCC) is the primary managing body of the
PypeIt project.  Its responsibilities include:

- Maintaining the health of the project, as assessed by the continued use of the
  code and timely attention to issues and pull requests.
    
- Recruitment for and management of the :ref:`pypeit_team`

- Establishing high-level development milestones

- Advertising opportunities for and coordinating project-related funding
  proposals and requests

The PCC has four members:

#. **Project Scientist (PS)**: The chair of the PCC, who guides the overall
   scientific direction of the project, elected for a term of three years.

#. **Developer Representative (DR)**: Presents and advocates for priorities of
   the developer community, elected for a term of two years.

#. **End-User Representative (UR)**: Presents and advocates for priorities of
   the end-user community, elected for a term of two years.

#. **Advisory Council Representative (AR)**: Presents and advocates for
   priorities of the broader astronomical community, both those that do and do
   not currently use PypeIt.  The Advisory Council Representative follows an
   independent selection process, as outlined by the
   :ref:`advisory-council-charter`.  While on the PCC, the PAC Representative is
   considered a Voting Member; however, they lose that designation once they
   vacate their spot on the PCC.

The PCC meets at least quarterly.  Decisions are made via consensus following
discussions that are recorded in openly accessible meeting notes.  If consensus
cannot be reached or when the PCC thinks it is otherwise appropriate or
necessary, decisions will be brought to a vote by the
:ref:`governance-voting-members`.

Nominations for PCC positions are held yearly, depending on the positions that
are reaching the end of their term.  Nominations can only come from
:ref:`governance-voting-members`, and Voting Members can self-nominate.  The
current PCC will vet the nominees and contact each nominee to ensure they accept
the nomination.  See further information :ref:`here<governance-appointments>`.
The procedure for removal of PCC members is discussed
:ref:`here<governance-removals>`.

Developer and End-User Representatives
++++++++++++++++++++++++++++++++++++++

Ultimately, the success of PypeIt will rely on **developers** who are willing
and able to improve and maintain the code and **end-users** that can
successfully use it to extract data for direct scientific analysis from raw
observations.  PypeIt encourages significant overlap between these two groups;
i.e., we encourage end-users to develop the code and developers to use PypeIt
for science projects.  The purpose of having Developer and End-User
Representatives on the PC, however, is to ensure that the needs specific to each
group are well represented within the project leadership.

Each Representative is expected to:

- Facilitate the collection of information from developers/end-users regarding
  issues with PypeIt, ranging from minor bugs to significant feature
  deficiencies

- As members of the PCC, present and advocate for development effort needed to
  address those needs

In addition, the DR will encourage development toward the milestones developed
by the PCC, and the UR will encourage community-driven efforts to test and
collect feedback on specific new/improved functionality that stems from each
milestone.

.. _governance-advisory:

PypeIt Advisory Council (PAC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PypeIt must continue to be a valuable resource for the astronomical community.
The PypeIt Advisory Council is a largely self-governed body whose primary task
is to ensure that the PypeIt Team is aware of the needs of the community and to
provide independent feedback on PypeIt's impact.  The PAC charter is kept
:ref:`here<advisory-council-charter>`, separate from this governance document.

.. _governance-maintain:

Maintainers
~~~~~~~~~~~

PypeIt Maintainers are a subset of developers that will have direct write access
to the PypeIt repository.  Maintainers are expected to:

- Provide user support, parsing feedback into GitHub Issues, as necessary

- Develop new features, perform code refactors, and address issues raised

- Improve test coverage

- Ensure the documentation is up-to-date

- Review pull requests, and potentially execute the suite of tests needed for
  pull-request reviews

Pull requests must be approved by at least one maintainer; however,
these maintainers may delegate the review to other community members without
write privileges.  In that case, the approval message from the maintainer (the
requirement for the merge to continue) must indicate that they are basing their
approval on the designated review.

Although maintainers are not designated as experts in specific parts of the
code, the PCC should make every effort to ensure that the group of maintainers,
as a whole, provides expertise spanning the full code base.

Maintainers are appointed by consensus of the PCC, or by vote if the PCC does
not come to consensus.  Maintainers are not appointed for a fixed term.
However, as a maintainer's priorities change, they may be asked (by the PCC or
by vote) to relinquish their write access if they are unable or unwilling to
continue fulfilling the duties of the role.  Maintainers can voluntarily
relinquish their write privileges at anytime.  The procedure for removal of a
Maintainer is discussed :ref:`here<governance-removals>`.

.. _governance-voting-members:

Voting Members
~~~~~~~~~~~~~~

Voting Members are members of the PypeIt community that have demonstrated their
commitment to the success of the project via significant contributions.  We
acknowledge contributions take many forms, including but not limited to
participating in discussions in our Users Slack Workspace; reporting issues to
our GitHub repository; submitting pull requests with small bug fixes,
documentation or testing improvements, or large feature improvements; and
participating in project maintenance and governance.  The PypeIt team should
maintain no fewer than 10 Voting Members.

Nominations for voting members can be sent to the PCC at any time for
consideration.  The PCC will vet the nominees and contact the nominee to ensure
they accept the nomination.  Self-nominations are acceptable.  Assuming both the
PCC and the nominee accept the nomination, voting member nominations will be
considered alongside PCC elections once per year.  Voting members are not
elected to a fixed term.  The procedure for removal of a Voting Members is
discussed :ref:`here<governance-removals>`.

.. _governance-other:

Other Named Roles
~~~~~~~~~~~~~~~~~

- **Ombudsperson**: The Ombudsperson provides a confidential point of contact
  for code-of-conduct violations or related concerns.  They also may help
  mediate disputes within the PypeIt team, as appropriate.  The Ombudsperson is
  nominated and appointed without a fixed term following the same procedure as
  for the PCC members.  The Ombudsperson cannot be a member of the PCC.

.. _governance-voting:

Voting Procedures
-----------------

For any vote to be valid, at least two thirds of the Voting Members must
participate.  Members that do not participate in any given vote are assumed to
have abstained from voting and will not be included in the vote tally.

Unless an off-cycle election is necessary, elections for PCC and Voting Members
are held yearly.  Voting Members are expected to vote in at least one of every
two yearly elections.  Voting members that have not participated in consecutive
elections will be asked to assume emeritus status, losing their active voting
privileges.  Participation in additional bespoke votes/elections do not count
for or against this tally.

Once a vote is called, votes are received over a two week period by anonymous
ballot.  Votes are tallied and reported by the PCC or the Ombudsperson.  Once
reported, the results of the vote take immediate effect, except for changes to
PCC membership (see below).

.. _governance-appointments:

Elections and Appointment Procedures
------------------------------------

Nominees for any named position are successfully elected if (1) they receive the
most votes when multiple nominees are under consideration or (2) they receive a
simple majority of votes when they are the sole nominee.

Appointments of any elected nominee take immediate effect, except for PCC
positions.  For PCC positions, there must be a transition period of at least one
month.

None of PypeIt's named positions are term limited; however, PCC members have
fixed terms where they must be re-elected to their positions.  The terms of each
PCC position is listed below.

=====================================  ======
Position                               Term
=====================================  ======
Project Scientist (PS)                 3 yr
Developer Representative (DR)          2 yr
End-User Representative (UR)           2 yr
Advisory Council Representative (AR)   ...
=====================================  ======

To maintain continuity, elections of the DR and UR are staggered.  The term of
the AR is set by the PAC but should be staggered with the PS.

If a member of the PCC is removed, their replacement will complete the original
term and must be re-elected to the position following the regular election
cycle.

.. _governance-removals:

Removal and Modification Procedures
-----------------------------------

Removal of individuals from named positions (except for Maintainers) must follow
a motion from a Voting Member that is seconded.  If the vote does not involve a
member of the PCC, the PCC will officiate the vote, as usual; otherwise, the
Ombudsperson will officiate.  Similarly, to modify this governance charter, the
modification must be proposed and seconded by Voting Members.

Removal of an individual from a position or modification of this charter
requires (1) at least two thirds of the Voting Members participate (as always)
and (2) at least three quarters of the tallied votes be in favor of the motion.

For Maintainers, removal can be proposed by a single PCC member.  All PCC
members must participate in the vote.  If the vote is split, the motion is taken
to the full list of Voting Members.

Institution
-----------

Instituting this new governance model will proceed as follows:

- The Ombudsperson through this process will be John O'Meara.

- The current PypeIt leadership team (Prochaska, Hennawi, Westfall, and Holden)
  will nominate potential Voting Members, Maintainers, and PCC members (except
  for the PAC Representative).  The current PypeIt leadership team members and
  the Ombudsperson will be de facto Voting Members.

- Voting Member nominees will be approached and either agree or disagree to
  accept the nomination.  Those that agree will be instituted as Voting Members
  immediately.

- Led by the current PypeIt leadership team, all Voting Members will iterate on
  and ultimately vote to ratify this governance model.

- Voting members will be able to make additional nominations for the first PCC,
  following the accepted procedure.  All PCC nominees will be approached, and if
  they accept the nomination, they will be part of the first election.  

- Once all the nominations are finalized, the first election will be held to
  fill the Project Scientist, Developer Representative, and End-User
  Representative positions on the PCC.

- The newly elected PCC will work with NOIRLab to setup the first PypeIt
  Advisory Council and iterate on the :ref:`advisory-council-charter`.

- The newly formed PAC will nominate their first representative for the PCC, and
  both the PAC Charter and the first PAC Representative will be voted on by all
  Voting Members.

- The now fully formed PCC will finalize the list of Maintainers, following the
  accepted procedure.

