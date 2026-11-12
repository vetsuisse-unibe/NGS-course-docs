-- Pandoc Lua filter for the NGS course exercises (HTML version 2).
-- Turns recurring patterns in the Markdown into styled callout boxes:
--   "Question:" / "Questions:" paragraph + following list  -> questions box
--   "#### To remember" + following list                     -> warning box
--   "**For the brave.** ..." paragraph + following code      -> challenge box
--   "#### Please ask for help ..."                           -> info box
-- It also drops the first level-1 heading (the page title is shown in the
-- header instead) and turns a final level-1 heading into a closing banner.

local stringify = pandoc.utils.stringify

local function title_div(text)
  return pandoc.Div({ pandoc.Plain({ pandoc.Str(text) }) }, pandoc.Attr("", { "callout-title" }))
end

local function callout(kind, title, blocks)
  local content = { title_div(title) }
  for _, b in ipairs(blocks) do table.insert(content, b) end
  return pandoc.Div(content, pandoc.Attr("", { "callout", "callout-" .. kind }))
end

local function is_list(b)
  return b and (b.t == "OrderedList" or b.t == "BulletList")
end

-- "Question:" or "Questions:" as the first word of a paragraph
local function question_para(b)
  if b.t ~= "Para" or #b.content == 0 then return false end
  local first = b.content[1]
  return first.t == "Str" and (first.text == "Question:" or first.text == "Questions:")
end

-- paragraph starting with bold "For the brave."
local function brave_para(b)
  if b.t ~= "Para" or #b.content == 0 then return false end
  local first = b.content[1]
  return first.t == "Strong" and stringify(first):match("^For the brave")
end

local function strip_label(para)
  local inl = para.content:clone()
  inl:remove(1)
  while #inl > 0 and (inl[1].t == "Space" or inl[1].t == "SoftBreak") do inl:remove(1) end
  if #inl == 0 then return nil end
  return pandoc.Para(inl)
end

function Pandoc(doc)
  local blocks = doc.blocks
  local out = pandoc.List()
  local first_h1_seen = false

  -- index of the last level-1 heading, if it is not the first one
  local last_h1 = nil
  for i, b in ipairs(blocks) do
    if b.t == "Header" and b.level == 1 then last_h1 = i end
  end

  local i = 1
  while i <= #blocks do
    local b = blocks[i]

    if b.t == "Header" and b.level == 1 and not first_h1_seen then
      first_h1_seen = true
      if i == last_h1 then last_h1 = nil end
      i = i + 1

    elseif b.t == "Header" and b.level == 1 and i == last_h1 then
      out:insert(pandoc.Div({ pandoc.Para(b.content) }, pandoc.Attr("", { "closing" })))
      i = i + 1

    elseif question_para(b) then
      local label = b.content[1].text == "Question:" and "Question" or "Questions"
      local inner = pandoc.List()
      local rest = strip_label(b)
      if rest then inner:insert(rest) end
      local j = i + 1
      while j <= #blocks do
        local nb = blocks[j]
        if is_list(nb) then
          inner:insert(nb); j = j + 1
        elseif nb.t == "CodeBlock" and (is_list(blocks[j + 1]) or #inner == 0) then
          -- code that belongs to a question (between list items)
          inner:insert(nb); j = j + 1
        else
          break
        end
      end
      out:insert(callout("questions", label, inner))
      i = j

    elseif b.t == "Header" and b.level >= 3 and stringify(b):match("^To remember") then
      local inner = pandoc.List()
      local j = i + 1
      while j <= #blocks and is_list(blocks[j]) do inner:insert(blocks[j]); j = j + 1 end
      out:insert(callout("remember", "Remember", inner))
      i = j

    elseif b.t == "Header" and b.level >= 3 and stringify(b):match("^Please ask") then
      out:insert(callout("info", "Help", { pandoc.Para(b.content) }))
      i = i + 1

    elseif brave_para(b) then
      local inner = pandoc.List()
      local rest = strip_label(b)
      if rest then inner:insert(rest) end
      local j = i + 1
      while j <= #blocks and blocks[j].t == "CodeBlock" do inner:insert(blocks[j]); j = j + 1 end
      out:insert(callout("brave", "For the brave", inner))
      i = j

    else
      out:insert(b)
      i = i + 1
    end
  end

  doc.blocks = out
  return doc
end
