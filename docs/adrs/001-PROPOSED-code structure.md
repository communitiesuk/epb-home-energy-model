# General Code Structure

|||
|-----|-----|
| Sequence | 001 |
| Date | 2026-02-13 |
| Author| Brian Jones |
| Authors Job Title | Lead Technical Architect |
| Approver | |
| Approvers job title ||
| Status | **PROPOSED** / ACCEPTED / REJECTED / SUPERSEDED / DEPRECIATED|

## Decision

Accept God Class anti-pattern for now.

## Context

This repository is a library which speeds up the calculation of energy costs as prescribed  in [a source repository](https://dev.azure.com/BreGroup/_git/Home%20Energy%20Model) maintained by [BRE Group](https://bregroup.com/)

The source repository is written by subject matter experts for subject matter expert understanding and communication. It is authoritative, but not performant.

This repository is written and maintained by software engineers, who are not in contact with the subject matter experts, The aim is to maximise performance with authority remaining with the slower source repository. 

This decision is a consequence of Conways Law `The structure of code mirrors the structure of the business`. The business in this context being two separate government departments. [Department for Energy Security & Net Zero](https://www.gov.uk/government/organisations/department-for-energy-security-and-net-zero) who maintain the calculation, and [Ministry Housing Communities & Local Government](https://www.gov.uk/government/organisations/ministry-of-housing-communities-local-government) who maintain this library.  

### Mitigation

In order to mitigate the business issues some limited refactoring will take place. But the priority will still be for ease of mapping of changes to the source repository over pure sftware considerations.

## Consequences

The trade offs of maintaining god classes are:
- Slower onboarding for new developers coming onto the team.
  - New engineers onboarding into the team will need significantly longer to absorb understanding of this repository.
  - When making changes that involve the classes affected the velocity of feature delivery will be reduced.
- Easier to map between original repository and this library
  - This repo is about performance, which it executes admirably, but the bridge of understanding between what is in the original.
  - The source repository is under development by subject matter experts, with frequent updates which are likely to continue in the long term. Maintaining the god classes in a similar manner to the source material reduces the 'mapping complexity' for developers.


## Alternatives Considered

### Applying Domain Driven Design Principles

Applying Domain Driven Design methodology may give a higher quality software product. However this becomes increasingly difficult without access to subject matter expertise to understand the domain completely. 

## Links

- [God class anti pattern reference](https://docs.embold.io/anti-patterns/#god-class)
- [Source Repository](https://dev.azure.com/BreGroup/_git/Home%20Energy%20Model)