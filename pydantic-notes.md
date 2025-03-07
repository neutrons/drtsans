# DRTSANS - Pydantic Notes

## Unrelated

- Recommend increasing minimum python version from 3.8 to 3.10 if not 3.12
- Typing:
  - Python 3.12 introduces the `type` keyword which would be handy for defining and reusing custom types in the configuration model
  - If this is too new, Python 3.10 introduced the `TypeAlias` type hint which would allow similar functionality
- Could be using Enums for things like centering method and instrument names
  - Would require some refactoring to handle things like going from enum to string and back,
    or cases where the entry in the config.json is lower case but the enum is upper case

## Linkml

- Would in theory be perfect for this use case, but possible it won't work because of the complex "multiple types" like `Union[int, List[int]]`, or `Union[str, int]`

  - There does exist a [`union_of`](https://linkml.io/linkml-model/latest/docs/union_of/) property in linkml, but it's unclear how it would handle these cases and interact with having to manage multiple type-specific constraints
  - linkml is more strict with typing and assumes a "range" is either a single type, a class, or an enum, which in theory would be better but may not be backwards compatible with some existing configuration files

- As an aside, I've contributed to the LinkML project, and they do have a Slack where we could ask questions if needed and work with the developers to address any issues that come up or features that would be beneficial.

## Pydantic

### Pros

- Includes [built-in validation](https://docs.pydantic.dev/latest/concepts/validators/) and error handling
  - Also allows for [custom validation](https://docs.pydantic.dev/latest/examples/custom_validators/)
    - See also: [Reusing custom validators](https://blog.det.life/pydantic-for-experts-reusing-importing-validators-2a4300bdcc8u)
- Lots of upfront work, but would increase maintainability and readability of the codebase in the long run, improving the developer experience
- Integrates well with `Annotated` and `Field` for more complex types, allowing for inclusion of metadata and default values
  - Pydantic has built-in tools to extract the metadata from an Annotated field and use it as validator(s).
- Allows for [custom data types](https://docs.pydantic.dev/latest/concepts/types/#custom-data-types) to be defined and reused
- By default, attempts to coerce values to the correct type, which would be useful for backwards compatibility with existing configuration files
  - Includes [Strict Mode](https://docs.pydantic.dev/latest/concepts/strict_mode/) for more control over this behavior
- As classes, attributes can be added to the models to include metadata, such as descriptions, default values, etc.
  - Can be easily accessed via dot notation and include type-hints for better readability
  - Can also be dumped to a dictionary/JSON for serialization

### Cons / Concerns / Areas of Effort

- Significant work to convert the existing json schema to pydantic models
  - Additionally, would need to refactor the existing code to use these models rather than the existing dict objects
- Significant manual testing to ensure that the pydantic models are working as expected and are backwards compatible with the existing configuration files
- Would need to write custom validators for some of the more complex types, such as strings that have to match a pattern, or numbers that have to be within a certain range
  - Though much of this could be done with the built-in validators or by reusing existing code in `drtsans.redparams`
- ReductionParameters would contain several nested models, which makes it difficult to create from a nested dictionary
  - Would need to write a custom `from_dict` method to handle this
- May need to write a custom `json` method to handle the serialization of the models
