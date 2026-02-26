# General Code Structure

|||
|-----|-----|
| Sequence | 001 |
| Date | 2026-02-13 |
| Author| Brian Jones |
| Authors Job Title | Lead Technical Architect |
| Approver | |
| Approver's job title ||
| Status | **PROPOSED** / ACCEPTED / REJECTED / SUPERSEDED / DEPRECIATED|

## Decision

Accept God Class anti-pattern for now.

## Context

This repository is a library that speeds up the calculation of energy costs as prescribed in [a source repository](https://dev.azure.com/BreGroup/_git/Home%20Energy%20Model) maintained by [BRE Group](https://bregroup.com/)

The source repository is written by subject matter experts for subject matter expert understanding and communication, as opposed to execution on a chip. It is authoritative, but not performant.

This repository is written and maintained by software engineers who are not in contact with the subject matter experts. The aim is to maximize performance with authority remaining with the slower source repository. 

This decision is a consequence of Conway's Law: `The structure of code mirrors the structure of the business`. The business in this context consists of two separate government departments. [Department for Energy Security & Net Zero](https://www.gov.uk/government/organisations/department-for-energy-security-and-net-zero), which maintains the calculation, and the [Ministry of Housing, Communities & Local Government](https://www.gov.uk/government/organisations/ministry-of-housing-communities-local-government), which maintains this library.  

### Mitigation

To mitigate the business issues, some limited refactoring will take place. However, the priority will remain on ease of mapping changes to the source repository over pure software considerations.

## Consequences

The trade-offs of maintaining god classes are:
- Slower onboarding for new developers joining the team.
  - New engineers onboarding to the team will need significantly longer to understand this repository.
  - When making changes involving these classes, the velocity of feature delivery will be reduced.
- Easier to map between the original repository and this library.
  - This repository prioritizes performance, which it achieves admirably, and maintaining structural similarity to the original creates a bridge of understanding.
  - The source repository is under active development by subject matter experts with frequent updates likely to continue long term. Maintaining the god classes in a similar manner to the source material reduces the 'mapping complexity' for developers.


## Alternatives Considered

### Applying Domain Driven Design Principles

Applying Domain-Driven Design methodology may produce a higher quality software product. However, this becomes increasingly difficult without access to subject matter expertise to fully understand the domain. This approach is unlikely to succeed without partnership. Partnership is unlikely because the goals of the two groups involved are not well aligned. Additionally, this would require a sustained relationship between the two entities over which we have no control.

## Links

- [God class anti pattern reference](https://docs.embold.io/anti-patterns/#god-class)
- [Source Repository](https://dev.azure.com/BreGroup/_git/Home%20Energy%20Model)