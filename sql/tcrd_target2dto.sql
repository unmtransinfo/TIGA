SELECT
	t.id target_id,
	p.id protein_id,
	p2dto.dtoid,
	p2dto.generation dto_generation,
	dto.name dto_class
FROM
	target t
	JOIN t2tc ON t2tc.target_id = t.id
	JOIN protein p ON p.id = t2tc.protein_id
	LEFT OUTER JOIN p2dto ON p2dto.protein_id = p.id
	JOIN dto ON dto.dtoid = p2dto.dtoid
ORDER BY
	t.id, p.id, p2dto.generation DESC
	;
