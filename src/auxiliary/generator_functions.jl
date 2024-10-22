# function point_data_field end

# function generate_point_data_field_functions(fields)
#     pdf_funcs = Expr[]
#     for field in fields
#         field_symb = QuoteNode(field)
#         pdf_func = quote
#             function Peridynamics.point_data_field(s::AbstractStorage, ::Val{$(field_symb)})
#                 return getfield(s, $(field_symb))
#             end
#         end
#         push!(pdf_funcs, pdf_func)
#     end
#     return Expr(:block, pdf_funcs...)
# # end

# macro point_data_field(storage, field)
#     if storage isa Symbol || !(storage isa Expr && storage.head === :.)
#         msg = "argument `$storage` is not a valid storage input!\n"
#         throw(ArgumentError(msg))
#     end
#     if !(field isa QuoteNode && field.value isa Symbol)
#         msg = "argument `$field` is not a valid field input!\n"
#         throw(ArgumentError(msg))
#     end
#     local _field_symb = field isa QuoteNode ? field : QuoteNode(field)
#     # if hasmethod(Peridynamics.point_data_field, Tuple{storage,Val(field)})
#     #     println("Hello! This method exists")
#     # else
#     #     println("-------")
#     # end
#     local _pdf_func = quote
#         function Peridynamics.point_data_field(storage::$(esc(storage)),
#                                                ::Val{$(_field_symb)})
#             return Peridynamics.getfield(storage, $(_field_symb))
#         end
#     end
#     return _pdf_func
# end
